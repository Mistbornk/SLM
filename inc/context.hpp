#pragma once
#include <string>
#include <iostream>
#include <string>
#include <unordered_map>

struct Options {
    std::string bam_path;
    std::string fasta_path;
    std::string output_path;
    size_t num_threads;
    size_t kmer;
    size_t bin_size;
};

enum class ReferenceGenome { GRCh37, GRCh38 };

struct ParRegion {
    size_t start;
    size_t end;
};

using ParMap = std::unordered_map<std::string, std::vector<ParRegion>>;

const std::unordered_map<ReferenceGenome, ParMap> PAR_REGIONS = {
    {
        ReferenceGenome::GRCh37, {
            {"X", {{60001, 2699520}, {154931044, 155260560}}},
            {"Y", {{10001, 2649520}, {59034050, 59363566}}},
            {"chrX", {{60001, 2699520}, {154931044, 155260560}}},
            {"chrY", {{10001, 2649520}, {59034050, 59363566}}}
        }
    },
    {
        ReferenceGenome::GRCh38, {
            {"X", {{10001, 2781479}, {155701383, 156030895}}},
            {"Y", {{10001, 2781479}, {56887903, 57217415}}},
            {"chrX", {{10001, 2781479}, {155701383, 156030895}}},
            {"chrY", {{10001, 2781479}, {56887903, 57217415}}}
        }
    }
};

constexpr const double PSEUDO_COUNT = 1.0;

enum class Sex { Male, Female};

inline bool verbose = false;
constexpr const char* INDEX_NAME = "index.bin";
constexpr const char* INTERVAL_NAME = "interval.bed";
constexpr const char* DISCARD_NAME = "discard.bed";
constexpr const char* READCOUNT_NAME = "read_count.bed";
constexpr const char* LOG2RATIO_NAME = "log2ratio.bed";

inline void verbose_log(const std::string& message) {
    if (verbose) {
        std::cout << message << std::endl;
    }
}