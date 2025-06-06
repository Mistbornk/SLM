#pragma once
#include <htslib/sam.h>
#include <string>
#include <vector>
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

std::vector<std::pair<std::string, int>> get_chromosomes_and_lengths(bam_hdr_t* header);
