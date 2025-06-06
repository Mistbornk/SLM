#pragma once
#include <biovoltron/algo/align/tailor/index.hpp>
#include <biovoltron/utility/istring.hpp>

#include "verbose.hpp"

#include <htslib/faidx.h>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <string_view>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <tuple>
#include <omp.h>

using namespace biovoltron;
using istring = std::basic_string<ichar>;
namespace fs = std::filesystem;

struct chr_info {
	std::string name;
	size_t len;
	size_t total_len;
};
using KmerResult = std::tuple<std::vector<chr_info>, std::vector<uint8_t>>;

KmerResult count_unique_kmers(const std::string& fasta_path, const std::string& index_path, int kmer_size, const std::string& output_path, int num_threads = 8);




