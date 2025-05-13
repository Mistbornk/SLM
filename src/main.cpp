#include "slmseg.hpp"
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
    SLMSeg<double> segtest;
    string file = "/Users/mistborn/Desktop/VScode/C++/SLM/testdata";
    string chr = "1";
    segtest.load_signal_file(file, chr);
    segtest.SLM();

    std::vector<std::vector<double>> data_bulk = segtest.data_seg();
    std::vector<double> data_seg = data_bulk[0];
    std::vector<unsigned int> pos_data = segtest.pos_data();
    
    std::ofstream segout("chr1_segment.txt");
    std::string chrom = "chr1";
    
    double current_value = data_seg[0];
    int start_index = 0;
    
    for (size_t i = 1; i < data_seg.size(); ++i) {
        if (data_seg[i] != current_value) {
            segout << chrom << "\t"
                   << pos_data[start_index] << "\t"
                   << (pos_data[i - 1]) << "\t"
                   << std::fixed << std::setprecision(6) << current_value << "\n";
            current_value = data_seg[i];
            start_index = i;
        }
    }
    segout << chrom << "\t"
           << pos_data[start_index] << "\t"
           << (pos_data.back()) << "\t"
           << std::fixed << std::setprecision(6) << current_value << "\n";
    
    segout.close();
    
    std::cout << "Save to csv successfully.\n";


    SLMSeg<double> segtest2;
    file = "/Users/mistborn/Desktop/VScode/C++/SLM/testdata";
    chr = "2";
    segtest2.load_signal_file(file, chr);
    segtest2.SLM();

    data_bulk = segtest2.data_seg();
    data_seg = data_bulk[0];
    pos_data = segtest2.pos_data();
    
    std::ofstream segout2("chr2_segment.txt");
    chrom = "chr2";
    
    current_value = data_seg[0];
    start_index = 0;
    
    for (size_t i = 1; i < data_seg.size(); ++i) {
        if (data_seg[i] != current_value) {
            segout2 << chrom << "\t"
                   << pos_data[start_index] << "\t"
                   << (pos_data[i - 1]) << "\t"
                   << std::fixed << std::setprecision(6) << current_value << "\n";
            current_value = data_seg[i];
            start_index = i;
        }
    }
    segout2 << chrom << "\t"
           << pos_data[start_index] << "\t"
           << (pos_data.back()) << "\t"
           << std::fixed << std::setprecision(6) << current_value << "\n";
    
    segout2.close();
    
    std::cout << "Save to csv successfully.\n";

    return 0;
}
