#include "slmseg.hpp"

void calSLM(std::string input_path, size_t num_threads) {
    const std::string out_path = "/home/max/SLM/SLM/data/all_segments.bed";
    std::vector<std::string> chromosomes;

    for (int i = 1; i <= 22; ++i) chromosomes.push_back(std::to_string(i));
    chromosomes.push_back("X");
    chromosomes.push_back("Y");

    size_t N = chromosomes.size();

    std::vector<std::vector<double>> all_seg(N);
    std::vector<std::vector<Bin_Region>> all_pos(N);

    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (size_t i = 0; i < N; ++i) {
        const std::string& chr = chromosomes[i];
        SLMSeg<double> hslm;
        hslm.load_signal_file(input_path, chr);
        hslm.HSLM();

        auto data_bulk = hslm.data_seg();
        auto pos_data = hslm.pos_data();  // 假設有此函式：bin 對應位置

        if (data_bulk.empty() || pos_data.empty()) continue;

        all_seg[i] = data_bulk[0];
        all_pos[i] = pos_data;

        verbose_log("Segmented chr" + chr);
    }
    std::ofstream segout(out_path);
    segout << std::fixed << std::setprecision(6);

    for (size_t i = 0; i < N; ++i) {
        const std::string& chrom = chromosomes[i];
        const auto& data_seg = all_seg[i];
        const auto& pos_data = all_pos[i];
        if (data_seg.empty() || pos_data.empty()) continue;

        double current_value = data_seg[0];
        int start_index = 0;

        for (size_t j = 1; j < data_seg.size(); ++j) {
            if (data_seg[j] != current_value) {
                segout << chrom << "\t"
                       << pos_data[start_index].start << "\t"
                       << pos_data[j - 1].start + pos_data[j - 1].len << "\t"
                       << current_value << "\n";
                current_value = data_seg[j];
                start_index = j;
            }
        }
        // 輸出最後一段
        segout << chrom << "\t"
               << pos_data[start_index].start << "\t"
               << pos_data.back().start + pos_data.back().len << "\t"
               << current_value << "\n";
    }

    segout.close();
}

