#pragma once
#include <vector>
#include <biovoltron/algo/align/tailor/index.hpp>
#include <biovoltron/utility/istring.hpp>

#include "context.hpp"
#include "fasta2kmer.hpp"

#include <htslib/faidx.h>
#include <algorithm>
#include <vector>
#include <tuple>
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

struct Intv {
    std::string name;
	// [start, end)
    size_t start;
    size_t end;
    double count = 0;
    //double gc_percent;
};
using BinningResult = std::tuple<std::vector<Intv>, std::vector<Intv>>;

BinningResult Binning_with_unique_kmer(const Options& option, const KmerResult& kmer_result);
std::vector<Intv> load_Bedfile_interval(const std::string& file_path);