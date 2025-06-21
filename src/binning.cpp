#include "binning.hpp"

std::vector<Intv> load_Bedfile_interval(
    const std::string& output_path
) {
    std::string file_path = output_path + INTERVAL_NAME;

    std::ifstream ifs(file_path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open interval input file: " << file_path << "\n";
        exit(1);
    }

    std::vector<Intv> bins_intvs;
    std::string line;

    while (std::getline(ifs, line)) {
        std::istringstream iss(line);
        std::string name;
        size_t start, end;

        if (!(iss >> name >> start >> end)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            exit(1);
        }

        Intv intv{name, start, end, 0};  // count 預設為 0

        bins_intvs.push_back(intv);
    }

    ifs.close();
    verbose_log("Successfully loaded intervals from " + file_path);

    return bins_intvs;
}

void write_Bedfile_interval(
	const std::string& output_path,
	const std::vector<Intv>& bins_intvs,
	const std::vector<Intv>& discard_intvs
) {
	// write bins intervals
	std::string interval_path = output_path + INTERVAL_NAME;
	std::ofstream ofs_interval(interval_path);
	if (!ofs_interval.is_open()) {
		std::cerr << "Failed to open interval output file: " << interval_path << "\n";
		exit(1);
	}

	for (const auto& intv : bins_intvs) {
		ofs_interval << intv.name << "\t" << intv.start << "\t" << intv.end << "\n";
	}

	ofs_interval.close();
	verbose_log("Successfully write intervals to " + interval_path);

	// write discard intervals
	std::string discard_path = output_path + DISCARD_NAME;
	std::ofstream ofs_discard(discard_path);
	if (!ofs_discard.is_open()) {
		std::cerr << "Failed to open discard output file: " << discard_path << "\n";
		exit(1);
	}

	for (const auto& intv : discard_intvs) {
		ofs_discard << intv.name << "\t" << intv.start << "\t" << intv.end << "\n";
	}

	ofs_discard.close();
	verbose_log("Successfully write intervals to " + discard_path);
}

BinningResult Binning_with_unique_kmer(
	const Options& option,
	const KmerResult& kmer_result
) {
	const std::string output_path(option.output_path);
	const size_t unique_number{option.bin_size};
	const size_t num_threads{option.num_threads};

	std::vector<std::vector<Intv>> bins_per_chr;
	std::vector<std::vector<Intv>> discard_per_chr;
	const auto& [chr_data, non_unique] = kmer_result;
	bins_per_chr.resize(chr_data.size());
	discard_per_chr.resize(chr_data.size());

	#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
	for (size_t k=0; k<chr_data.size(); ++k) {

		size_t start{chr_data[k].total_len};
		size_t end{start + chr_data[k].len};
		std::string name(chr_data[k].name);

		size_t count = 0;
		size_t total_len{chr_data[k].total_len};

		std::vector<Intv> temp_bins;
		std::vector<Intv> temp_discards;

		for (size_t i = start, j = start; j < end; ++j) {
			if (!non_unique[j]) ++count;
			
			if (count >= unique_number) {
				if (j - i + 1 <= 2 * unique_number) {
					temp_bins.emplace_back(name, i - total_len, j - total_len);
				}
				i = j;
				count = 0;
			}
			if (j - i + 1 > 2 * unique_number) {
				temp_discards.emplace_back(name, i - total_len, j - total_len);
				i = j;
				count = 0;
			}
		}
		bins_per_chr[k] = std::move(temp_bins);
		discard_per_chr[k] = std::move(temp_discards);
	}

	// flatten
	std::vector<Intv> bins_intvs_flat;
	std::vector<Intv> discard_intvs_flat;

	for (size_t k = 0; k < chr_data.size(); ++k) {
		bins_intvs_flat.insert(bins_intvs_flat.end(), bins_per_chr[k].begin(), bins_per_chr[k].end());
		discard_intvs_flat.insert(discard_intvs_flat.end(), discard_per_chr[k].begin(), discard_per_chr[k].end());
	}

	verbose_log("Binning Success.");
	// write_Bedfile_interval(output_path, bins_intvs_flat, discard_intvs_flat);

	return std::make_tuple(bins_intvs_flat, discard_intvs_flat);
}
