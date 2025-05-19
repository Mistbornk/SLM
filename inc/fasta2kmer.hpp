#ifndef FASTA2KMER_H
#define FASTA2KMER_H
#include <htslib/faidx.h>
#include <unordered_map>
#include <algorithm>

class kmer_info {
public:
	std::string chr;
	unsigned int pos;
};

void count_unique_kmers(const std::string& fasta_path, int kmer_size, const std::string& output_path, int num_threads = 8);

#endif



