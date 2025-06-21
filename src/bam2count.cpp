#include "bam2count.hpp"
#include <atomic>
std::atomic<bool> printed_first_zero(false);

void BAM_INFO::bam_file_init()
{
	// open bam
	bam_file = sam_open(file_path.c_str(), "r");
    if (!bam_file) {
        std::cerr << "Can not open bam file\n";
        std::exit(1);
    }

	// read header
	header = sam_hdr_read(bam_file);
    if (!header) {
        std::cerr << "Can not read header\n";
        sam_close(bam_file);
        std::exit(1);
    }

	// load index
	idx = sam_index_load(bam_file, file_path.c_str());
    if (!idx) {
        std::cerr << "Can't find BAM index (use samtools to construct .bai)\n";
        bam_hdr_destroy(header);
        sam_close(bam_file);
        std::exit(1);
    }
}

void BAM_INFO::bam_file_close()
{
    sam_close(bam_file);
    bam_hdr_destroy(header);
	hts_idx_destroy(idx);
}

static double median_of_vector(std::vector<double> v) {
    size_t n = v.size();
    if (!n) throw std::runtime_error("Cannot median empty vector");
    
    size_t mid = n / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    double m = v[mid];
    if (!(n&1)) {
        // even: avg with lower median
        std::nth_element(v.begin(), v.begin() + mid - 1, v.begin() + mid);
        m = (m + v[mid - 1]) * 0.5;
    }
    return m;
}

inline bool is_in_par(const std::string& chr, size_t start, size_t end, ReferenceGenome ref) {
    auto it_genome = PAR_REGIONS.find(ref);
    if (it_genome == PAR_REGIONS.end()) return false;

    auto it_chr = it_genome->second.find(chr);
    if (it_chr == it_genome->second.end()) return false;

    for (const auto& region : it_chr->second) {
        if (!(end < region.start || start > region.end)) {
            return true;  // overlap
        }
    }
    return false;
}

void print_par_regions(
    const std::vector<std::vector<Intv>>& bins_by_chr,
    const std::vector<std::string>& chr_names,
    ReferenceGenome ref_version
) {
    for (size_t i = 0; i < bins_by_chr.size(); ++i) {
        const auto& chr = chr_names[i];
        const auto& bin_list = bins_by_chr[i];

        for (const auto& intv : bin_list) {
            if (is_in_par(chr, intv.start, intv.end, ref_version)) {
                std::cout << chr << "\t"
                          << intv.start << "\t"
                          << intv.end << "\t"
                          << std::fixed << std::setprecision(4)
                          << intv.count << "\n";
            }
        }
    }
}

double get_baseline_multiplier(const std::string& chr, Sex sex) {
    if (chr == "X" || chr == "chrX") return (sex == Sex::Male) ? 0.5 : 1.0;
    if (chr == "Y" || chr == "chrY") return (sex == Sex::Male) ? 0.5 : 0.0; // female 無 Y
    return 1.0;
}

inline bool is_autosome(const std::string& chr) {
    std::string name = chr;

    if (name.rfind("chr", 0) == 0)
        name = name.substr(3);

    try {
        int n = std::stoi(name);
        return (n >= 1 && n <= 22);
    } catch (...) {
        return false;
    }
}

void normalization(
    std::vector<std::vector<Intv>>& bins_by_chr,
    std::vector<std::string>& chr_names,
    const size_t num_threads,
    Sex sample_sex,
    ReferenceGenome sample_ref
) {
    size_t N = bins_by_chr.size();
    std::vector<double> contig_medians_over_total;
    contig_medians_over_total.reserve(N);
    
    for (size_t i = 0; i < N; ++i) {
        if (!is_autosome(chr_names[i])) continue; 

        const auto& bin_list = bins_by_chr[i];
        std::vector<double> raw_count;
        raw_count.reserve(bin_list.size());

        std::transform(
            bin_list.begin(),
            bin_list.end(),
            std::back_inserter(raw_count),
            [](const auto& intv) { return intv.count; }
        );
        
        double med = median_of_vector(raw_count);
        double total = std::accumulate(raw_count.begin(), raw_count.end(), 0.0);
        if (total > 0.0) {
            contig_medians_over_total.emplace_back(med / total);
        } else {
            std::cerr << "[WARN] total count for chromosome " << chr_names[i] << " is zero. Skipping.\n";
        }
    }

    double baseline = median_of_vector(contig_medians_over_total);

    std::vector<double> log2_ratios;
    for (size_t i = 0; i < N; ++i) {
        const auto& bin_list = bins_by_chr[i];
        const auto& chr = chr_names[i];

        double multiplier = get_baseline_multiplier(chr, sample_sex);

        for (const auto& intv : bin_list) {
            double local_multiplier = multiplier;
            if (is_in_par(chr, intv.start, intv.end, sample_ref)) {
                local_multiplier = 1.0;
            }

            if (local_multiplier == 0.0) {
                log2_ratios.push_back(NAN);
            } else {
                double copy_ratio = (intv.count + PSEUDO_COUNT) / ((baseline * local_multiplier));
                double log2r = std::log2(copy_ratio);
                log2_ratios.push_back(log2r);
            }
        }
    }

    std::vector<double> valid_log2;
    valid_log2.reserve(log2_ratios.size());

    std::copy_if(
        log2_ratios.begin(),
        log2_ratios.end(),
        std::back_inserter(valid_log2),
        [](double v) { return !std::isnan(v); }
    );

    double median_log2 = median_of_vector(valid_log2);

    // Write back centered log2(copy-ratio)
    size_t idx = 0;
    for (auto& bin_list : bins_by_chr) {
        for (auto& intv : bin_list) {
            double v = log2_ratios[idx++];
            intv.count = std::isnan(v) ? NAN : (v - median_log2);
        }
    }
}

void gc_bias_corrections(
    const std::string& fasta_path,
    std::vector<Intv>& bins_intvs,
    const size_t num_threads,
    int gc_bin_size_pct = 1
) {
    size_t N = bins_intvs.size();
    std::vector<double> raw_counts, gc_percents;
    raw_counts.resize(N);
    gc_percents.resize(N);

    #pragma omp parallel num_threads(num_threads)
    {
        // load local fatsa for every threads
        faidx_t* local_fai = fai_load(fasta_path.c_str());
        if (!local_fai) {
            #pragma omp critical
            std::cerr << "Failed to open fai index in thread.\n";
            exit(1);
        }

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < N; ++i) {
            raw_counts[i] = bins_intvs[i].count;

            // load  interval seq
            hts_pos_t len = 0;
            char* seq = faidx_fetch_seq64(local_fai, bins_intvs[i].name.c_str(), 
                                        bins_intvs[i].start, bins_intvs[i].end - 1, &len);
            if (!seq || size_t(len) != bins_intvs[i].end - bins_intvs[i].start) {
                if (seq) free(seq);
                fai_destroy(local_fai);
                std::cerr << "Error: failed to fetch sequence for " << bins_intvs[i].name 
                        << ":" << std::to_string(bins_intvs[i].start) << "-" 
                        << std::to_string(bins_intvs[i].end) << "\n";
                continue;
            }

            // count gc
            size_t gc = 0;
            for (size_t j = 0; j < size_t(len); ++j) {
                char c = seq[j];
                if (c=='G' || c=='C' || c=='g' || c=='c') ++gc;
            }
            free(seq);

            gc_percents[i] = double(gc) / double(len) * 100.0;
            //bins_intvs[i].gc_percent = gc_percents[i];
        }

        fai_destroy(local_fai);
    }

    // put interval count into gc percent slots
    int num_slots = 100 / gc_bin_size_pct + 1;
    std::vector<std::vector<double>> bin2counts;
    bin2counts.resize(num_slots);
    for (size_t i = 0; i < N; ++i) {
        int bin = std::min(int(gc_percents[i] / gc_bin_size_pct), num_slots - 1);
        bin2counts[bin].emplace_back(raw_counts[i]);
    }

    // cal every gc slots median
    std::vector<double> slot_medians(num_slots, 0.0);
    for (int s = 0; s < num_slots; ++s) {
        auto& v = bin2counts[s];
        if (!v.empty()) slot_medians[s] = median_of_vector(v);
    }

    // cal all intervals' global median
    double global_med = median_of_vector(raw_counts);

    // write back to Intc.count
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (size_t i = 0; i < N; ++i) {
        int bin = std::min(int(gc_percents[i] / gc_bin_size_pct), num_slots - 1);
        double m = slot_medians[bin];

        // interpolation
        if (m == 0.0) {
            int left = bin - 1;
            int right = bin + 1;  

            while (left >= 0 && slot_medians[left] == 0.0) --left;
            while (right < num_slots && slot_medians[right] == 0.0) ++right;

            if (left >= 0 && right < num_slots) {
                double x1 = left * gc_bin_size_pct;
                double x2 = right * gc_bin_size_pct;
                double y1 = slot_medians[left];
                double y2 = slot_medians[right];
                double x = gc_percents[i];

                m = ((x2 != x1)) ? y1 + (y2 - y1) * (x - x1) / (x2 - x1) : (y1 + y2) / 2.0;

            } else if (left >= 0) {
                m = slot_medians[left];
            } else if (right < num_slots) {
                m = slot_medians[right];
            } else {
                // 完全找不到 → fallback 預設值
                m = global_med;
            }
        }
        // correct the raw counts by gc
        bins_intvs[i].count = raw_counts[i] / m * global_med;;
    }
}

std::vector<size_t> transfer_tid2chr(
    BAM_INFO& bam_info,
    const std::vector<CHR_INFO>& chr_data
) {
    std::unordered_map<std::string, size_t> name2chr;
    name2chr.reserve(chr_data.size());

    for (size_t i = 0; i < chr_data.size(); ++i) {
        name2chr[chr_data[i].name] = i;
    }

    int n_target = bam_info.header->n_targets;
    std::vector<size_t> tid2chr(n_target, std::string::npos);

    for (int tid = 0; tid < n_target; ++tid) {
        std::string_view bam_chr(bam_info.header->target_name[tid]);
        auto it = name2chr.find(std::string(bam_chr));
        if (it == name2chr.end()) {
            throw std::runtime_error("Chromosome in BAM header not found in FASTA: " + std::string(bam_chr));
        }
        tid2chr[tid] = it->second;
    }

    # if 0
    for (int tid = 0; tid < n_target; ++tid) {
        size_t chr_idx = tid2chr[tid];
        if (chr_idx == std::string::npos) {
            std::cout << "tid " << tid
                    << " (" << bam_info.header->target_name[tid] << ")"
                    << " -> NOT FOUND in chr_data\n";
        } else {
            std::cout << "tid " << tid
                    << " (" << bam_info.header->target_name[tid] << ")"
                    << " -> chr_data[" << chr_idx << "] = "
                    << chr_data[chr_idx].name << "\n";
        }
    }
    # endif

    return tid2chr;
}

//void write_Bedfile_ReadCount(
//	const std::string& output_path,
//	const std::vector<Intv>& bins_intvs
//) {
//	std::string readcount_path = output_path + READCOUNT_NAME;
//	std::ofstream ofs_readcount(readcount_path);
//	if (!ofs_readcount.is_open()) {
//		std::cerr << "Failed to open interval output file: " << readcount_path << "\n";
//		exit(1);
//	}

//    for (const auto& intv : bins_intvs) {
//        # if 0
//        ofs_readcount << intv.name << "\t" 
//                      << intv.start << "\t" 
//                      << intv.end << "\t" 
//                      << std::fixed << std::setprecision(4)
//                      << intv.count << "\t" 
//                      << intv.gc_percent << "\n";
//        # endif

//        # if 1
//        ofs_readcount << intv.name << "\t" 
//                      << intv.start << "\t" 
//                      << intv.end << "\t"
//                      << intv.count << "\n";
//        # endif 
//    }

//	ofs_readcount.close();
//	verbose_log("Successfully write bins read count to " + readcount_path);
//}

void write_Bedfile_ReadCount(
	const std::string& output_path,
	const std::vector<std::vector<Intv>>& bins_intvs
) {
	std::string readcount_path = output_path + LOG2RATIO_NAME;
	std::ofstream ofs_readcount(readcount_path);
	if (!ofs_readcount.is_open()) {
		std::cerr << "Failed to open interval output file: " << readcount_path << "\n";
		exit(1);
	}
    for (const auto& bin_list : bins_intvs) {
        for (const auto& intv : bin_list) {
            ofs_readcount << intv.name << "\t" 
                        << intv.start << "\t" 
                        << intv.end << "\t"
                        << intv.count << "\n";
        }
    }
	ofs_readcount.close();
	verbose_log("Successfully write bins read count to " + readcount_path);
}

inline bool is_proper(const bam1_t* aln) 
{
    // p.22
    // “target counts” are summarized as the number of reads falling in such intervals that are:  
    // primary alignments, not duplicates, properly paired, in forward orientation, with MAPQ ≥ 3
    return !(aln->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP)) &&
           (aln->core.flag & BAM_FPROPER_PAIR) &&
           !(aln->core.flag & BAM_FREVERSE) &&
           aln->core.qual >= 3;
}

int count_reads_in_interval(
    const BAM_INFO& bam_info, 
    const Intv& intv,
    const std::vector<size_t>& tid2chr,
    const std::vector<CHR_INFO>& chr_data, 
    const std::vector<uint8_t>& non_unique
) { 

    // name to tid
    int tid = bam_name2id(bam_info.header, intv.name.c_str());
    if (tid < 0) {
        throw std::runtime_error("Chromosome name not found: " + intv.name);
    }

	// query interval
    hts_itr_t* iter = sam_itr_queryi(bam_info.idx, tid, intv.start, intv.end);
    if (!iter) {
        throw std::runtime_error(
            "Failed to query BAM for " + intv.name + ":" +
            std::to_string(intv.start) + "-" + std::to_string(intv.end)
        );
    }

    // count reads
    bam1_t* aln = bam_init1();
    int count = 0;
    const size_t chr_idx = tid2chr[tid];
    const size_t base_off = chr_data[chr_idx].total_len; 
    while (sam_itr_next(bam_info.bam_file, iter, aln) >= 0) {

        if (!is_proper(aln)) continue; 
        // p.22 starting on a k-mer unique position
        size_t pos = aln->core.pos;
        size_t global_index = base_off + pos;
        if ( non_unique[global_index]|| global_index >= non_unique.size()) continue;

        ++count;
    }

    // close and destroy
    bam_destroy1(aln);
    hts_itr_destroy(iter);

    return count;
}

std::vector<std::vector<Intv>> count_reads_in_bam(
    const Options& option, 
    std::vector<Intv>& bins_intvs, 
    const KmerResult& kmer_result,
    Sex sample_sex,
    ReferenceGenome sample_ref
) {
    const std::string output_path(option.output_path);
    const std::string bam_file_path(option.bam_path);
    const std::string fasta_path(option.fasta_path);
    const size_t num_threads{option.num_threads};

    const auto& [chr_data, non_unique] = kmer_result;

    BAM_INFO bam_info(bam_file_path);
    const auto tid2chr = transfer_tid2chr(bam_info, chr_data);

    verbose_log("Counting bins reads ...");
    #pragma omp parallel num_threads(num_threads)
    {
        BAM_INFO bam_info(bam_file_path);

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < bins_intvs.size(); ++i) {
            auto& intv = bins_intvs[i];
            intv.count = count_reads_in_interval(bam_info, intv, tid2chr, chr_data, non_unique);
        }
    }
    verbose_log("Count Reads Success.");

    verbose_log("Correcting GC content ...");
    gc_bias_corrections(fasta_path, bins_intvs, num_threads, 1);

    #if 0
        write_Bedfile_ReadCount(output_path, bins_intvs);
    #endif

    auto is_valid_chrom = [](const std::string& chr) {
        std::string c = chr;
        if (c.rfind("chr", 0) == 0) c = c.substr(3);

        if (c == "X" || c == "Y") return true;

        try {
            int n = std::stoi(c);
            return (n >= 1 && n <= 22);
        } catch (...) {
            return false;
        }
    };

    std::vector<std::vector<Intv>> bins_by_chr;
    std::vector<std::string> chr_names;
    std::vector<Intv> current_chr_group;
    std::string current_chr = bins_intvs[0].name;
    
    for (const auto& intv : bins_intvs) {
        if (intv.name != current_chr) {
            if (is_valid_chrom(current_chr)) {
                if (current_chr == "Y" || current_chr == "chrY") {
                    current_chr_group.erase(
                        std::remove_if(
                            current_chr_group.begin(),
                            current_chr_group.end(),
                            [&](const Intv& i) {
                                return is_in_par(i.name, i.start, i.end, sample_ref);
                            }
                        ),
                        current_chr_group.end()
                    );
                }
                bins_by_chr.push_back(std::move(current_chr_group));
                chr_names.push_back(current_chr);
            }
            current_chr = intv.name;
            current_chr_group.clear();
        }
        current_chr_group.push_back(intv);
    }
    if (is_valid_chrom(current_chr)) {
        if (current_chr == "Y" || current_chr == "chrY") {
            current_chr_group.erase(
                std::remove_if(
                    current_chr_group.begin(),
                    current_chr_group.end(),
                    [&](const Intv& i) {
                        return is_in_par(i.name, i.start, i.end, sample_ref);
                    }
                ),
                current_chr_group.end()
            );
        }
        bins_by_chr.push_back(std::move(current_chr_group));
        chr_names.push_back(current_chr);
    }

    verbose_log("normalized...");
    normalization(bins_by_chr, chr_names, num_threads, sample_sex, sample_ref);

    write_Bedfile_ReadCount(output_path, bins_by_chr);

    return bins_by_chr;
}