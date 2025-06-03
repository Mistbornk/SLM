#include "fasta2kmer.hpp"

std::vector<chr_info> chr_data;
//std::unordered_map<std::string, unsigned int> chrtoid;
std::vector<uint8_t> non_unique; 
istring concat_iref;
//std::vector<std::uint64_t> kmers;
const std::string index_name("index.bin");

void save_non_unique(const std::string& output_path, const std::vector<uint8_t>& non_unique) {
    // 壓縮成 bit vector
    std::vector<uint8_t> packed((non_unique.size() + 7) / 8, 0);
    for (size_t i = 0; i < non_unique.size(); ++i) {
        if (non_unique[i]) {
            packed[i / 8] |= (1 << (i % 8));
        }
    }

    // 確保路徑存在（選擇性）
    std::filesystem::create_directories(std::filesystem::path(output_path).parent_path());

    // 打開檔案寫 binary
    std::ofstream ofs(output_path, std::ios::binary);
    if (!ofs) {
        std::cerr << "Failed to open file for writing: " << output_path << "\n";
        exit(1);
    }
    ofs.write(reinterpret_cast<const char*>(packed.data()), packed.size());
    ofs.close();

    std::cout << "Saved non-unique bitvector to " << output_path << "\n";
}

void make_FMindex(const std::string& file_path) 
{   
    std::cout << "  -Making FM index ...\n";

    auto index = FMIndex<>{}; 
    index.build(concat_iref);
    std::cout << "  -Index construct success\n";

    std::string out_path = file_path + index_name;
    std::ofstream out(out_path, std::ios::binary); 
    index.save(out);
    out.close();
    std::cout << "  -Save success\n";
}


void count_unique_kmers(const std::string& fasta_path, const std::string& index_path, 
    int kmer_size, const std::string& output_path, int num_threads
) {

	faidx_t* fai = fai_load(fasta_path.c_str());
	if (!fai) {
        std::cerr << "Failed to load fasta index." << std::endl;
        exit(1);
    }

    //istring ref = Codec::to_istring("ATCCGTCGGAATCCGTCGGAATCCGTCGGAATCCGTCGGAATCCGTCGGA");
    //auto index = FMIndex<>{.LOOKUP_LEN = 5}; 
    //index.build(ref);
    //std::cerr << "Success\n";

    //std::string s = "ATCCGTCGGA";
    //auto read = Codec::to_istring(s);
    //auto [begin, end, offset] = index.get_range(read, 0);

    //std::cerr << "Read: " << s << "\n";
    //std::cerr << "Begin: " << begin << ", End: " << end << ", Offset: " << offset << "\n";
    //std::cerr << "Count: " << end - begin << "\n";
    //std::cerr << "=====================" << "\n";
    //if (end - begin > 1) {
    //    auto hits = index.get_offsets(begin, end);
    //    for (auto& hit : hits) {
    //        std::cout << "Pos: " << hit << "\n";
    //    }
    //}

    // read chrom name and len 
    int chr_num = faidx_nseq(fai);
    for (size_t  i = 0, total_len = 0; i < chr_num; ++i) {
        std::string chr_name(faidx_iseq(fai, i));
        if (chr_name.empty()) {
            std::cerr << "Warning: empty chromosome name at index " << i << "\n";
            continue;
        }

        int chr_len = faidx_seq_len(fai, chr_name.c_str());
        if (chr_len <= 0) {
            std::cerr << "Warning: failed to get length for " << chr_name << "\n";
            continue;
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
    for (int i = 0 ; i < chr_num; ++i) {

        std::string& chr_name = chr_data[i].name;
        size_t chr_len = chr_data[i].len;
        hts_pos_t len = 0;
        size_t total_len = chr_data[i].total_len;
        
        char* seq = faidx_fetch_seq64(fai, chr_name.c_str(), 0, chr_len - 1, &len);
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
        free(seq);
        concat_iref.append(iseq);

        std::cout << "Read " << chr_name << "\n";
    }
    fai_destroy(fai);

    // make and save FM index
    if (!fs::exists(index_path + index_name)) {
        std::cerr << "Need to construct FM index first ...\n"; 
        make_FMindex(index_path);
    }

    auto index = FMIndex<>{};
    {
        auto ifs = std::ifstream{index_path + index_name};
        index.load(ifs);
    }

    std::cout << "Start to Proccess Chromosome\n";
    std::cout << "Kmer size: " << kmer_size << "\n";
    size_t chr_idx = 0; 
    size_t chr_end = chr_data[chr_idx].total_len + chr_data[chr_idx].len;
    //std::cout << "start " << chr_data[chr_idx].name << " " << chr_end << "\n";
    //for (auto& info : chr_data) {
    //    std::cout << info.name << ", " << info.total_len << ", " << info.len <<"\n";
    //}
    //std:: cout << concat_iref.size() << "\n";

    for (size_t j = 0; j <= concat_iref.size() - kmer_size; ++j) {   
        if (non_unique[j]) continue;

        istring_view ikmer_view(concat_iref.data() + j, kmer_size);
        const auto [begin, end, offset] = index.get_range(ikmer_view, 0);
        size_t count = end - begin;

        if (count > 1) {
            auto hits = index.get_offsets(begin, end);
            for (auto& hit : hits) {
                if (hit < non_unique.size()) {
                    non_unique[hit] = 1;
                }
            }
        }

        if (j >= chr_end) {
            std::cout << "Processed " << chr_data[chr_idx].name << "\n";
            ++chr_idx;
            if (chr_idx < chr_data.size()) {
                chr_end = chr_data[chr_idx].total_len + chr_data[chr_idx].len;
            }
        }
    }

    std::cout << "Successfully done. Now output the non unique vector\n";
    save_non_unique(output_path, non_unique);

}


std::string get_chr_name_from_pos(size_t pos) {
    for (size_t i = 0; i < chr_data.size(); ++i) {
        if (i+1 < chr_data.size()) {
            if (pos >= chr_data[i].total_len && pos < chr_data[i+1].total_len) {
                return chr_data[i].name;
            }
        } else {
            if (pos >= chr_data[i].total_len) {
                return chr_data[i].name;
            }
        }
    }
    return "unknown";
}
