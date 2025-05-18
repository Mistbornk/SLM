#ifndef BAM2COUNT_H
#define BAM2COUNT_H
#include <string>
#include <htslib/sam.h>
#include "bam_utils.hpp"

int count_reads_in_interval(BAM_INFO* bam_info, const std::string& region_str);

#endif
