#pragma once
#include <string>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

struct Options {
    std::string bam_path;
    std::string fasta_path;
    std::string output_path;
    size_t num_threads;
    size_t kmer;
    size_t bin_size;
};

inline bool verbose = false;
constexpr const char* INDEX_NAME = "index.bin";
constexpr const char* INTERVAL_NAME = "interval.bed";
constexpr const char* DISCARD_NAME = "discard.bed";
constexpr const char* READCOUNT_NAME = "read_count.bed";

inline void verbose_log(const std::string& message) {
    if (verbose) {
        std::cout << message << std::endl;
    }
}