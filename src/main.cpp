#include <biovoltron/algo/align/tailor/index.hpp>
#include <biovoltron/utility/istring.hpp>

#include "slmseg.hpp"
#include "bam2count.hpp"
#include "bam_utils.hpp"
#include "fasta2kmer.hpp"
#include "binning.hpp"
#include "verbose.hpp"

#include <htslib/sam.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>


#define Kmer 35

int main(int argc, char* argv[]) {
    std::string file_path = "/mnt/users/philip/workshop/_share/20250217_dragen/hs37d5_graph/HG002.novaseq.pcr-free.30x/HG002.novaseq.pcr-free.30x.bam";  // ← 替換成你的路徑
    std::string fa_path = "/home/max/SLM/SLM/data/hs37d5.fa";
    std::string index_path = "/home/max/SLM/SLM/data/";
    std:: string output_path = "/home/max/SLM/SLM/data/";

    // thread and verbose
    int num_threads = std::thread::hardware_concurrency();
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--verbose") {
            verbose = true;
        } else {
            try {
                num_threads = std::stoi(arg);
            } catch (const std::invalid_argument&) {
                std::cerr << "Invalid argument: " << arg << std::endl;
                return 1;
            }
        }
    }

    auto kmer_result = count_unique_kmers(fa_path, index_path, Kmer, output_path, num_threads);

    Binning_with_unique_kmer(output_path, kmer_result, 1000);

    return 0;
}
