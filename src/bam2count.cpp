#include "bam2count.hpp"
#include <iostream>
#include <stdexcept>

inline bool is_primary_unique(const bam1_t* aln) {
    // p.22
    // “target counts” are summarized as the number of reads falling in such intervals that are: primary alignments, not duplicates, properly paired, in forward orientation, with MAPQ ≥ 3
    return !(aln->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP)) &&
           aln->core.qual >= 3;
}

std::mutex cout_mutex;

int count_reads_in_interval(BAM_INFO* bam_info, const std::string& region_str) 
{
	// query interval
	hts_itr_t* iter = sam_itr_querys(bam_info->idx, bam_info->header, region_str.c_str());
    if (!iter) {
        std::cerr << "Fail to query " << region_str << "\n";
        std::exit(-1);
    }

    // count reads
    bam1_t* aln = bam_init1();
    int count = 0;
    while (sam_itr_next(bam_info->bam_file, iter, aln) >= 0) {
        if (is_primary_unique(aln)) count++;
    }

    // close
    bam_destroy1(aln);
    hts_itr_destroy(iter);

    return count;
}

void thread_work(const std::string& chr, int len, const std::string& file_path) {
    BAM_INFO bam_info(file_path);  // 每個 thread 自己開 BAM

    std::string region = chr + ":1-" + std::to_string(len);
    int count = count_reads_in_interval(&bam_info, region);

    std::lock_guard<std::mutex> lock(cout_mutex);
    if (count >= 0)
        std::cout << chr << "\tlength=" << len << "\treads=" << count << "\n";
    else
        std::cerr << "read failed on " << region << "\n";
}

void count_reads_in_bam(const std::string& file_path, int num_threads) {
    BAM_INFO bam_info(file_path);  // 每個 thread 自己開 BAM

    auto chromosome_info = get_chromosomes_and_lengths(bam_info.header);

    // Create a vector of futures to hold the results
    std::vector<std::future<void>> futures;

    // Launch threads to count reads for each chromosome
    for (const auto& [chr, len] : chromosome_info) {
        futures.push_back(std::async(std::launch::async, thread_work, chr, len, file_path));
    }

    // Wait for all threads finish
    for (auto& future : futures) {
        future.get();
    }
}