#include "bam2count.hpp"
#include <iostream>
#include <stdexcept>

int count_reads_in_interval(BAM_INFO* bam_info, const std::string& region_str) 
{
	// query interval
	hts_itr_t* iter = sam_itr_querys(bam_info->idx, bam_info->header, region_str.c_str());
    if (!iter) {
        std::cerr << "Fail to query " << region_str << "\n";
        delete bam_info;
        std::exit(-1);
    }

    // count reads
    bam1_t* aln = bam_init1();
    int count = 0;
    while (sam_itr_next(bam_info->bam_file, iter, aln) >= 0) {
        if (!(aln->core.flag & BAM_FUNMAP))  // 排除未比對 reads
            count++;
    }

    // close
    bam_destroy1(aln);
    hts_itr_destroy(iter);

    return count;
}