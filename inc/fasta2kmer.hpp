#pragma once
#include <biovoltron/algo/align/tailor/index.hpp>
#include <biovoltron/utility/istring.hpp>

#include "context.hpp"

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

struct CHR_INFO {
	std::string name;
	size_t len;
	size_t total_len;
};
using KmerResult = std::tuple<std::vector<CHR_INFO>, std::vector<uint8_t>>;

KmerResult count_unique_kmers(const Options& option);




