#pragma once
#include <string>
#include <htslib/sam.h>
#include "bam_utils.hpp"
#include <future>
#include <thread>
#include <mutex>

void count_reads_in_bam(const std::string& file_path, int num_threads);

