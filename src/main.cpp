#include <biovoltron/algo/align/tailor/index.hpp>
#include <biovoltron/utility/istring.hpp>

#include "slmseg.hpp"
#include "bam2count.hpp"
#include "fasta2kmer.hpp"
#include "binning.hpp"
#include "context.hpp"

#include <htslib/sam.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>

Options parse_args(int argc, char* argv[]) 
{
    Options option;
    option.num_threads = std::thread::hardware_concurrency();
    option.kmer = 35;
    option.bin_size = 1000;

    for (int i = 1; i < argc; ++i) {
        std::string key(argv[i]);

        try {
            if (key == "-verbose") {
                verbose = true;
            }
            else if (key == "-fasta" || key == "-f") {
                if (i + 1 >= argc) throw std::invalid_argument("Missing value");
                option.fasta_path = argv[++i];
            }
            else if (key == "-bam" || key == "-b") {
                if (i + 1 >= argc) throw std::invalid_argument("Missing value");
                option.bam_path = argv[++i];
            } 
            else if (key == "-output" || key == "-o") {
                if (i + 1 >= argc) throw std::invalid_argument("Missing value");
                std::string output_path(argv[++i]);

                if (output_path.empty()) {
                    throw std::invalid_argument("Output path cannot be empty.");
                }

                option.output_path = (output_path.back() == '/') ? output_path : output_path + '/';
            }
            else if (key == "-thread" || key == "-t") {
                if (i + 1 >= argc) throw std::invalid_argument("Missing value");
                option.num_threads = std::stoi(argv[++i]);
            }
            else if (key == "-kmer" || key == "-k") {
                if (i + 1 >= argc) throw std::invalid_argument("Missing value");
                option.kmer = std::stoi(argv[++i]);
            }
            else if (key == "-bin") {
                if (i + 1 >= argc) throw std::invalid_argument("Missing value");
                option.bin_size = std::stoi(argv[++i]);
            }
            else {
                std::cerr << "Unknown option: " << key << std::endl;
                exit(1);
            }
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "Invalid value for option: " << key << ", expected a value." << std::endl;
            exit(1);
        }
        catch (const std::out_of_range& e) {
            std::cerr << "Value out of range for option: " << key << std::endl;
            exit(1);
        }
        catch (...) {
            std::cerr << "Unknown error occurred while parsing option: " << key << std::endl;
            exit(1);
        }
    }

    if (option.fasta_path.empty() || option.output_path.empty()) {
        std::cerr << "Missing required options: -fasta, -bam, -output are required.\n";
        exit(1);
    }

    return option;
}

int main(int argc, char* argv[]) {
    Options option = parse_args(argc, argv);

    if (!fs::exists(option.output_path + INTERVAL_NAME)) {
        auto kmer_result = count_unique_kmers(option);
        auto [bins_intv, discard_intv] = Binning_with_unique_kmer(option, kmer_result);
        count_reads_in_bam(option, bins_intv, kmer_result);
    }
    else {
        auto kmer_result = count_unique_kmers(option);
        auto bins_intv = load_Bedfile_interval(option.output_path);
        count_reads_in_bam(option, bins_intv, kmer_result);
    }

    return 0;
}
