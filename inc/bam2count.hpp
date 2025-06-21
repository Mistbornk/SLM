#pragma once
#include <htslib/sam.h>

#include "binning.hpp"
#include "fasta2kmer.hpp"
#include "context.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <future>
#include <thread>
#include <mutex>
#include <utility>  // for std::pair

struct BAM_INFO {
	BAM_INFO (const std::string& bam_path) : file_path(bam_path) {
		bam_file_init();
	}

	~BAM_INFO() { 
		bam_file_close();
	}

    void bam_file_init();
    void bam_file_close();	

	samFile* bam_file;
	bam_hdr_t* header;
	hts_idx_t* idx;
	const std::string file_path;
};

std::vector<std::vector<Intv>> count_reads_in_bam(
	const Options& option, 
	std::vector<Intv>& bins_intvs, 
	const KmerResult& kemr_result,
    Sex sample_sex,
    ReferenceGenome sample_ref
);

