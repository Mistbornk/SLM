#include "slmseg.hpp"
#include "bam2count.hpp"
#include "bam_utils.hpp"
#include <htslib/sam.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <mutex>

int main() {
    std::string file_path = "/mnt/users/philip/workshop/_share/20250217_dragen/hs37d5_graph/HG002.novaseq.pcr-free.30x/HG002.novaseq.pcr-free.30x.bam";  // ← 替換成你的路徑
    std::string region = "1:10000-11000";                   // 測試區間
    
    BAM_INFO* bam_info = new BAM_INFO(file_path);

    auto chromosome_info = get_chromosomes_and_lengths(bam_info->header);

    for (const auto& [chr, len] : chromosome_info) {
        std::string region = chr + ":1-" + std::to_string(len);
        int count = count_reads_in_interval(bam_info, region);
        if (count >= 0)
            std::cout << chr << "\tlength=" << len << "\treads=" << count << "\n";
        else
            std::cerr << "讀取失敗\n";
    }

    delete(bam_info);
    return 0;
}
