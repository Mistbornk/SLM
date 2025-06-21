#include "slmseg.hpp"

void calSLM(std::string input_path, size_t num_threads) {
    const std::string out_path = "/home/max/SLM/SLM/data/all_segments.bed";
    std::vector<std::string> chromosomes;

    for (int i = 1; i <= 22; ++i) chromosomes.push_back(std::to_string(i));
    chromosomes.push_back("X");
    chromosomes.push_back("Y");

    size_t N = chromosomes.size();

    std::vector<std::vector<double>> data_seg(N);
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (size_t i = 0; i < N; ++i) {
        const std::string& chr = chromosomes[i];
        SLMSeg<double> hslm;
        hslm.load_signal_file(input_path, chr);
        hslm.HSLM();

        auto data_bulk = hslm.data_seg();
        if (data_bulk.empty()) continue;

        data_seg[i] = data_bulk[0];
        verbose_log("Segmented chr"+chr);
    }

    

}

