#include "binning.hpp"
#include <fstream>

void write_Bedfile_interval(
	const std::string& output_path,
	std::vector<chr_info>& chr_data,
	std::vector<std::vector<Intv>>& bins_intvs,
	std::vector<std::vector<Intv>>& discard_intvs
) {
	// 1. write bins intervals
	std::string interval_name = output_path + "/interval.bed";
	std::ofstream ofs_interval(interval_name);
	if (!ofs_interval.is_open()) {
		std::cerr << "Failed to open interval output file: " << interval_name << "\n";
		exit(1);
	}

	for (size_t i = 0; i < bins_intvs.size(); ++i) {
		const auto& chr = chr_data[i];
		for (const auto& intv : bins_intvs[i]) {
			ofs_interval << chr.name << "\t" << intv.start - chr.total_len << "\t" << intv.end - chr.total_len << "\n";
		}
	}
	ofs_interval.close();
	verbose_log("Successfully write intervals to " + interval_name);

	// 2. write discard intervals
	std::string discard_name = output_path + "/discard.bed";
	std::ofstream ofs_discard(discard_name);
	if (!ofs_discard.is_open()) {
		std::cerr << "Failed to open discard output file: " << discard_name << "\n";
		exit(1);
	}

	for (size_t i = 0; i < discard_intvs.size(); ++i) {
		const auto& chr = chr_data[i];
		for (const auto& intv : discard_intvs[i]) {
			ofs_discard << chr.name << "\t" << intv.start - chr.total_len << "\t" << intv.end - chr.total_len << "\n";
		}
	}
	ofs_discard.close();
	verbose_log("Successfully write intervals to " + discard_name);
}

void Binning_with_unique_kmer(
	const std::string& output_path, 
	KmerResult& kmer_result, 
	size_t unique_number
){
	std::vector<std::vector<Intv>> bins_intvs;
	std::vector<std::vector<Intv>> discard_intvs;
	auto& [chr_data, non_unique] = kmer_result;
	bins_intvs.reserve(chr_data.size());

	for (auto& chr : chr_data) {
		size_t start{chr.total_len};
		size_t end{start + chr.len};
		size_t count = 0;
		std::vector<Intv> temp1_intervals;
		std::vector<Intv> temp2_intervals;

		for (size_t i = start, j = start; j < end; ++j) {
			if (!non_unique[j]) ++count;
			
			if (count >= unique_number) {
				if (j - i + 1 <= 2 * unique_number) {
					temp1_intervals.emplace_back(i, j + 1);
				}
				i = j + 1;
				count = 0;
			}
			if (j - i + 1 > 2 * unique_number) {
				temp2_intervals.emplace_back(i, j + 1);
				i = j + 1;
				count = 0;
			}
		}
		bins_intvs.emplace_back(temp1_intervals);
		discard_intvs.emplace_back(temp2_intervals);
	}
	write_Bedfile_interval(output_path, chr_data, bins_intvs, discard_intvs);
}
