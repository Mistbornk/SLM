#include "fasta2kmer.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <future>
#include <mutex>

using KmerMap = std::unordered_map<uint64_t, unsigned int>;
using KmerLocMap = std::unordered_map<uint64_t, kmer_info>;

// k-mer in unsigned 128-bit
inline bool encode_kmer(const std::string& seq, int start, int k, uint64_t& result) {
    result = 0;
    for (int i = 0; i < k; ++i) {
        char c = std::toupper(seq[start + i]);
        result <<= 2;
        switch (c) {
            case 'A': result |= 0; break;
            case 'C': result |= 1; break;
            case 'G': result |= 2; break;
            case 'T': result |= 3; break;
            default: return false; // N base
        }
    }
    return true;
}

std::pair<KmerMap, KmerLocMap> process_chromosome(const std::string& fasta_path, const std::string& chr, int kmer_size) {
    faidx_t* fai = fai_load(fasta_path.c_str());
    if (!fai) {
        throw std::runtime_error("Failed to load fasta index");
    }

    int len;
    char* seq = faidx_fetch_seq(fai, chr.c_str(), 0, faidx_seq_len(fai, chr.c_str()) - 1, &len);
    std::string str(seq);
    free(seq);
    fai_destroy(fai);

    KmerMap counts;
    KmerLocMap locs;

    for (int j = 0; j <= len - kmer_size; ++j) {
        uint64_t kmer;
        if (!encode_kmer(str, j, kmer_size, kmer)) continue;

        if (counts[kmer] == 0) {
            locs[kmer] = {chr, static_cast<unsigned int>(j)};
        }
        counts[kmer]++;
    }

    std::cout << "Processed chromosome: " << chr << std::endl;
    return {counts, locs};
}

void count_unique_kmers(const std::string& fasta_path, int kmer_size, const std::string& output_path, int num_threads) 
{

    if (kmer_size > 32) {
        std::cerr << "Error: k-mer size exceeds 32" << std::endl;
        exit(1);
    }

	faidx_t* fai = fai_load(fasta_path.c_str());
	if (!fai) {
        std::cerr << "Failed to load fasta index." << std::endl;
        exit(1);
    }
	
    int chr_num = faidx_nseq(fai);
    std::vector<std::string> chr_names;
    for (int i = 0; i < chr_num; ++i) {
        chr_names.emplace_back(faidx_iseq(fai, i));
    }
    fai_destroy(fai);

    std::vector<std::future<std::pair<KmerMap, KmerLocMap>>> futures;
    for (const auto& chr : chr_names) {
        futures.push_back(std::async(std::launch::async, process_chromosome, fasta_path, chr, kmer_size));
    }

    KmerMap merged_counts;
    KmerLocMap merged_locs;
    std::mutex merge_mutex;
    for (auto& fut : futures) {
        auto [counts, locs] = fut.get();
        std::lock_guard<std::mutex> lock(merge_mutex);  // thread-safe

        for (const auto& [kmer, count] : counts) {
            if (merged_counts[kmer] == 0) {
                merged_locs[kmer] = locs.at(kmer);
            }
            merged_counts[kmer] += count;
        }

        counts.clear();
        locs.clear();
    }
    //std::ofstream out(output_path);
    //if (!out) {
    //    std::cerr << "Unable to create output file: " << output_path << std::endl;
    //    exit(1);
    //}

    //for (const auto& [kmer, count] : merged_counts) {
    //    if (count == 1) {
    //        const auto& [chr, pos] = merged_locs[kmer];
    //        out << chr << '\t' << pos << '\t' << (pos + kmer_size) << '\n';
    //    }
    //}
    //out.close();
}