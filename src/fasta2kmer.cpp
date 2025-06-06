#include "fasta2kmer.hpp"

void make_FMindex(const std::string& file_path, const istring& concat_iref) 
{   
    verbose_log("  -Making FM index ...");

    auto index = FMIndex<>{}; 
    index.build(concat_iref);
    verbose_log("  -Index construct success");

    std::string out_path = file_path;
    std::ofstream out(out_path, std::ios::binary); 
    index.save(out);
    out.close();

    verbose_log("  -Save success");
}


KmerResult count_unique_kmers(const std::string& fasta_path, const std::string& index_path, 
    int kmer_size, const std::string& output_path, int num_threads
) {

    std::vector<chr_info> chr_data;
    std::vector<uint8_t> non_unique; 
    istring concat_iref;
    const std::string index_name("index.bin");
    const std::string bitvector_name("non_unique.bin");

	faidx_t* fai = fai_load(fasta_path.c_str());
	if (!fai) {
        std::cerr << "Failed to load fasta index." << std::endl;
        exit(1);
    }

    // read chrom name and len 
    int chr_num = faidx_nseq(fai);
    for (size_t  i = 0, total_len = 0; i < chr_num; ++i) {
        std::string chr_name(faidx_iseq(fai, i));
        if (chr_name.empty()) {
            std::cerr << "Empty chromosome name at index " << i << "\n";
            exit(1);
        }

        int chr_len = faidx_seq_len(fai, chr_name.c_str());
        if (chr_len <= 0) {
            std::cerr << "Failed to get length for " << chr_name << "\n";
            exit(1);
        }
        chr_data.emplace_back(chr_name, chr_len, total_len);
        total_len += chr_len;
    }

    // construct a bit vector
    if (!chr_data.empty()) {
        auto& last = chr_data.back();
        non_unique.resize(last.total_len + last.len, 0);
    }else {
        std::cerr << "Failed to construct bit vector." << std::endl;
        exit(1);
    }

    // read sequence info form fasta
    std::vector<istring> temp_sequences(chr_num);
    #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
    for (int i = 0 ; i < chr_num; ++i) {

        std::string& chr_name = chr_data[i].name;
        size_t chr_len = chr_data[i].len;
        hts_pos_t len = 0;
        size_t total_len = chr_data[i].total_len;

        faidx_t* local_fai = fai_load(fasta_path.c_str());
        if (!local_fai) {
            #pragma omp critical
            std::cerr << "Failed to open fai index in thread.\n";
            exit(1);
        }

        char* seq = faidx_fetch_seq64(local_fai, chr_name.c_str(), 0, chr_len - 1, &len);
        if (!seq) {
            std::cerr << "Error: failed to fetch sequence for " << chr_name << "\n";
            continue;
        }else if (len != chr_len) {
            std::cerr << "Warning: fetched length mismatch for " << chr_name 
              << " (expected " << chr_len << ", got " << len << ")\n";
        }  

        istring iseq;
        {
            iseq.reserve(chr_len);
            for (size_t j = 0; j < chr_len; ++j) {
                ichar val = Codec::to_int(seq[j]);
                if (val == 4) {
                    val = 0;  // N -> A
                    non_unique[total_len + j] = 1;
                }
                iseq.push_back(val);
            }
        }
        temp_sequences[i] = std::move(iseq);
        free(seq);
        fai_destroy(local_fai);

        #pragma omp critical
        verbose_log("Read " + chr_name);
    }

    // 依照原順序 append
    for (int i = 0; i < chr_num; ++i) {
        concat_iref.append(temp_sequences[i]);
    }

    // make and save FM index
    if (!fs::exists(index_path + index_name)) {
       verbose_log("Need to construct FM index first ...");
        make_FMindex(index_path + index_name, concat_iref);
    }

    // load
    auto index = FMIndex<>{};
    {
        verbose_log("Loading FM index ...");
        auto ifs = std::ifstream{index_path + index_name};
        index.load(ifs);
    }

    verbose_log("Start to cal unique kmer ...");
    verbose_log("Kmer size: " + std::to_string(kmer_size));

    // search kmer unique
    for (auto& chr : chr_data) {
        if (chr.total_len - kmer_size > 0 && chr.total_len + kmer_size < non_unique.size()) {
            size_t start = chr.total_len - kmer_size + 1;
            size_t end = chr.total_len;
            for (size_t i = start; i<end; ++i) {
                non_unique[i] = 1;
            }
        }
    }

    #pragma omp parallel for schedule(dynamic, 1000)
    for (size_t j = 0; j <= concat_iref.size() - kmer_size; ++j) {   
        if (non_unique[j]) continue;

        istring_view ikmer_view(concat_iref.data() + j, kmer_size);
        const auto [begin, end, offset] = index.get_range(ikmer_view, 0);
        size_t count = end - begin;

        if (count > 1) {
            auto hits = index.get_offsets(begin, end);
            for (auto& hit : hits) {
                if (hit < non_unique.size()) {
                    //#pragma omp atomic write
                    non_unique[hit] = 1;
                }
            }
        }
    }

    verbose_log("Successfully done unique kmer search !");

    return std::make_tuple(chr_data, non_unique);
}
