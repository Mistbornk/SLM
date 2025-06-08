#include "bam_utils.hpp"
#include <iostream>
#include <stdexcept>

void BAM_INFO::bam_file_init()
{
	// open bam
	bam_file = sam_open(file_path.c_str(), "r");
    if (!bam_file) {
        std::cerr << "Can not open bam file\n";
        std::exit(1);
    }

	// read header
	header = sam_hdr_read(bam_file);
    if (!header) {
        std::cerr << "Can not read header\n";
        sam_close(bam_file);
        std::exit(1);
    }

	// load index
	idx = sam_index_load(bam_file, file_path.c_str());
    if (!idx) {
        std::cerr << "Can't find BAM index (use samtools to construct .bai)\n";
        bam_hdr_destroy(header);
        sam_close(bam_file);
        std::exit(1);
    }
}

void BAM_INFO::bam_file_close()
{
    sam_close(bam_file);
    bam_hdr_destroy(header);
	hts_idx_destroy(idx);
}

std::vector<std::pair<std::string, int>> get_chromosomes_and_lengths(bam_hdr_t* header) 
{
    std::vector<std::pair<std::string, int>> result;

    for (int i = 0; i < header->n_targets; ++i) {
        std::string chr_name = header->target_name[i];
        int chr_len = header->target_len[i];
        result.emplace_back(chr_name, chr_len);
    }

    return result;
}