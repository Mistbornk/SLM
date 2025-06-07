#pragma once
#include <vector>
#include <biovoltron/algo/align/tailor/index.hpp>
#include <biovoltron/utility/istring.hpp>

#include "verbose.hpp"
#include "fasta2kmer.hpp"

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
using KmerResult = std::tuple<std::vector<chr_info>, std::vector<uint8_t>>;

struct Intv {
	// [start, end)
    size_t start;
    size_t end;
};

void Binning_with_unique_kmer(const std::string& output_path, KmerResult& kmer_result, size_t kmer_number);
