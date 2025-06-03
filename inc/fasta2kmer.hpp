#pragma once
#include <biovoltron/algo/align/tailor/index.hpp>
#include <biovoltron/utility/istring.hpp>
#include <htslib/faidx.h>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <string_view>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <future>
#include <mutex>
using namespace biovoltron;
using istring = std::basic_string<ichar>;
namespace fs = std::filesystem;

//struct kmer_info {
//    uint64_t kmer;
//    int chr;
//    int pos;
//};

struct chr_info {
	std::string name;
	size_t len;
	size_t total_len;
};

void count_unique_kmers(const std::string& fasta_path, const std::string& index_path, int kmer_size, const std::string& output_path, int num_threads = 8);




