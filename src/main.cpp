#include "slmseg.hpp"
#include "bam2count.hpp"
#include "bam_utils.hpp"
#include "fasta2kmer.hpp"
#include <htslib/sam.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#define Kmer 30

int main(int argc, char* argv[]) {
    std::string file_path = "/mnt/users/philip/workshop/_share/20250217_dragen/hs37d5_graph/HG002.novaseq.pcr-free.30x/HG002.novaseq.pcr-free.30x.bam";  // ← 替換成你的路徑
    std::string fa_path = "/home/max/SLM/SLM/data/hs37d5.fa";

    BAM_INFO bam_info(file_path);
    auto chromosome_info = get_chromosomes_and_lengths(bam_info.header);

    // thread
    int num_threads = std::thread::hardware_concurrency();
    if (argc > 1) num_threads = std::stoi(argv[1]);

    //count_reads_in_bam(file_path, num_threads);

    count_unique_kmers(fa_path, Kmer, "/home/max/SLM/SLM/data/uniqued_kmer.txt");

    return 0;
}
