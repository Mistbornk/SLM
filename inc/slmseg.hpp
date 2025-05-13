#ifndef SLMSEG_H_
#define SLMSEG_H_
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numbers>
#include <set>
#include <stdexcept>

// TFloat is a floating point value - float or double
template<class TFloat>
class SLMSeg {
public:
    SLMSeg(TFloat omega = 0.1, TFloat eta = 0.0001, unsigned int stepeta = 1000000, unsigned int fw = 1, unsigned int mw = 1)
    : omega_(omega), eta_(eta), stepeta_(stepeta), fw_(fw), mw_(mw) { }
    ~SLMSeg() { };

    bool load_signal_file(const std::string& filename, const std::string& chromosome);

    std::vector<unsigned int> get_breaks(std::vector<unsigned int>& path);

    TFloat elnsum(TFloat a, TFloat b) {
        if (a > b)
            return a + std::log1p(std::exp(b - a));
        else
            return b + std::log1p(std::exp(a - b));
    }

    void transemisi_HSLM(const std::vector<TFloat>& etavec, int NCov, int K, int NExp, int T,
                    std::vector<TFloat>& G, std::vector<std::vector<TFloat>>& P, std::vector<std::vector<TFloat>>& Emission);

    void transemisi_SLM(const TFloat eta, int K, int NExp, int T,
                    std::vector<TFloat>& G, std::vector<std::vector<TFloat>>& P, std::vector<std::vector<TFloat>>& Emission);                  
                    
    void bioviterbii_HSLM(const std::vector<TFloat>& etav, const std::vector<std::vector<TFloat>>& P, 
                    const std::vector<std::vector<TFloat>>& Emission, int  T, int  K, std::vector<unsigned int>& path, 
                    std::vector<std::vector<unsigned int>>& psi);
    
    void bioviterbii_SLM(const std::vector<TFloat>& etav, const std::vector<std::vector<TFloat>>& P, 
                    const std::vector<std::vector<TFloat>>& Emission, int  T, int  K, std::vector<unsigned int>& path, 
                    std::vector<std::vector<unsigned int>>& psi); 
    
    void param_est_seq();
    void muk_est();
    void joint_seg();
    void joint_seg_in();
    void filter_seg();
    void seg_results();

    void SLM();
    void HSLM();

    void write_breakpoints_to_file(const std::string& filename) const {
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            std::cerr << "Unable to open file for writing: " << filename << "\n";
            return;
        }
    
        for (unsigned int bp : total_pred_break_filtered_) {
            outfile << bp << "\n";
        }
    
        outfile.close();
        std::cout << "Breakpoints written to " << filename << "\n";
    }

    std::vector<TFloat>  mi(){ return mi_; }
    std::vector<TFloat>  smu(){ return smu_; }
    std::vector<TFloat>  sepsilon(){ return sepsilon_; }
    std::vector<std::vector<TFloat>> muk(){ return muk_; }
    std::vector<unsigned int> total_pred_break(){ return total_pred_break_;}
    std::vector<unsigned int> total_pred_break_filtered(){ return total_pred_break_filtered_; }
    std::vector<unsigned int> pos_data(){ return pos_data_; }
    std::vector<std::vector<TFloat>> data_seg(){ return data_seg_; }

private:
    std::vector<TFloat> signal_data_;
    std::vector<unsigned int> pos_data_;
    std::vector<std::vector<TFloat>> data_matrix;

    std::vector<TFloat> mi_;
    std::vector<TFloat> smu_;
    std::vector<TFloat> sepsilon_;
    TFloat omega_;
    TFloat eta_;
    TFloat stepeta_;
    TFloat fw_;
    TFloat mw_;

    std::vector<unsigned int> total_pred_break_;
    std::vector<unsigned int> total_pred_break_filtered_;
    std::vector<std::vector<TFloat>> muk_;
    std::vector<std::vector<TFloat>> data_seg_;	
};

template<class TFloat>
std::vector<unsigned int> SLMSeg<TFloat>::get_breaks(std::vector<unsigned int>& path) {
    std::vector<unsigned int> breakpoints;

    if (path.empty())
        return breakpoints;

    int T = path.size();
    unsigned int last_state = path[0];
    breakpoints.push_back(0);  // 起點一定是斷點

    for (int t = 1; t < T; ++t) {
        if (path[t] != last_state) {
            breakpoints.push_back(t);  // 狀態切換時記錄斷點
            last_state = path[t];
        }
    }
    // 確保最後一段被處理
    if (breakpoints.back() != T)
    breakpoints.push_back(T);

    return breakpoints;
}       

template<class TFloat>
bool SLMSeg<TFloat>::load_signal_file(const std::string& file_name, const std::string& chromosome) 
{
    std::ifstream infile(file_name);
    if (!infile.is_open()) {
        std::cerr << "Unable to open file: " << file_name << "\n";
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream ss(line);
        std::string col0, col1, col2;
        if (!std::getline(ss, col0, '\t')) continue;
        if (!std::getline(ss, col1, '\t')) continue;
        if (!std::getline(ss, col2, '\t')) continue;

        try {
            std::string chr  = col0;
            if (chr == chromosome) {
                pos_data_.push_back(static_cast<unsigned int>(std::stoul(col1)));
                signal_data_.push_back(static_cast<TFloat>(std::stod(col2)));
            }
        } catch (...) {
            std::cerr << "Unexpected error " << "\n";
            continue;
        }
    }

    infile.close();
    return true;
}

template<class TFloat>
void SLMSeg<TFloat>::transemisi_HSLM(
    const std::vector<TFloat>& etavec,
    int NCov,
    int K,
    int NExp,
    int T,
    std::vector<TFloat>& G, 
    std::vector<std::vector<TFloat>>& P, 
    std::vector<std::vector<TFloat>>& Emission) 
{
    TFloat PI = std::numbers::pi;
    // Step 1: GVECT (a simplify gaussion log-likelihood)
    std::vector<TFloat> gvect(K, 0.0);
    for (int i = 0; i < K; ++i) {
        TFloat gsum = 0.0;
        for (int j = 0; j < NExp; ++j) {
            TFloat diff = muk_[j][i] - mi_[j];
            gsum += - (diff * diff) / (2 * smu_[j] * smu_[j]);
        }
        gvect[i] = gsum;
    }

    // Step 2: normalize GVECT (likelihood) to G (probability)
    TFloat norm = gvect[0];
    for (int i = 1; i < K; ++i)
        norm = elnsum(norm, gvect[i]);
    for (int i = 0; i < K; ++i)
        G[i] = gvect[i] - norm;

    // Step 3: transition matrix P [K][K][NCov]
    for (int cov = 0; cov < NCov; ++cov) {
        for (int i = 0; i < K; ++i) {
            for (int j = 0; j < K; ++j) {
                if (i == j)
                    P[i][j + cov * K] = elnsum(std::log(1 - etavec[cov]), std::log(etavec[cov]) + G[j]);
                else
                    P[i][j + cov * K] = std::log(etavec[cov]) + G[j];
            }
        }
    }

    // Step 4: emission log-likelihood [K x T]
    for (int k = 0; k < K; ++k) {
        for (int t = 0; t < T; ++t) {
            TFloat val = 0.0;
            for (int j = 0; j < NExp; ++j) {
                TFloat diff = data_matrix[j][t] - muk_[j][k];
                val += -0.5 * (diff * diff) / (sepsilon_[j] * sepsilon_[j]) -
                    std::log(std::sqrt(2.0 * PI) * sepsilon_[j]);
            }
            Emission[k][t] = val;
        }
    }  
}

template<class TFloat>
void SLMSeg<TFloat>::transemisi_SLM(
    const TFloat eta,
    int K,
    int NExp,
    int T,
    std::vector<TFloat>& G, 
    std::vector<std::vector<TFloat>>& P, 
    std::vector<std::vector<TFloat>>& Emission) 
{
    TFloat PI = std::numbers::pi;

    // Step 1: GVECT
    std::vector<TFloat> gvect(K, 0.0);
    for (int i = 0; i < K; ++i) {
        TFloat gsum = 0.0;
        for (int j = 0; j < NExp; ++j) {
            TFloat diff = muk_[j][i] - mi_[j];
            gsum += - (diff * diff) / (2 * smu_[j] * smu_[j]);
        }
        gvect[i] = gsum;
    }

    // Step 2: Normalize GVECT to G
    TFloat norm = gvect[0];
    for (int i = 1; i < K; ++i)
        norm = elnsum(norm, gvect[i]);
    for (int i = 0; i < K; ++i)
        G[i] = gvect[i] - norm;

    // Step 3: Transition matrix P (K x K)
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            if (i == j)
                P[i][j] = elnsum(std::log(1.0 - eta), std::log(eta) + G[j]);
            else
                P[i][j] = std::log(eta) + G[j];
        }
    }

    // Step 4: Emission likelihoods (K x T)
    for (int k = 0; k < K; ++k) {
        for (int t = 0; t < T; ++t) {
            TFloat val = 0.0;
            for (int j = 0; j < NExp; ++j) {
                TFloat diff = data_matrix[j][t] - muk_[j][k];
                val += -0.5 * (diff * diff) / (sepsilon_[j] * sepsilon_[j]) -
                       std::log(std::sqrt(2.0 * PI) * sepsilon_[j]);
            }
            Emission[k][t] = val;
        }
    }
}

template<class TFloat>
void SLMSeg<TFloat>::bioviterbii_HSLM(
    const std::vector<TFloat>& etav, 
    const std::vector<std::vector<TFloat>>& P, 
    const std::vector<std::vector<TFloat>>& Emission,
    int T, 
    int K, 
    std::vector<unsigned int>& path, 
    std::vector<std::vector<unsigned int>>& psi) 
{
    std::vector<std::vector<TFloat>> delta(K, std::vector<TFloat>(T, 0.0));
    
    // initialization
    for (int i = 0; i < K; ++i) {
        delta[i][0] = etav[i] + Emission[i][0];
        psi[i][0] = 0;
    }

    // recursion
    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < K; ++j) {
            TFloat nummax = delta[0][t - 1] + P[0][j + K * (t - 1)]; // P[0][j][t]
            unsigned int ind = 0;

            for (int i = 1; i < K; ++i) {
                TFloat score = delta[i][t - 1] + P[i][j + K * (t - 1)]; // P[i][j][t-1]
                if (score > nummax) {
                    nummax = score;
                    ind = i;
                }
            }

            delta[j][t] = nummax + Emission[j][t];
            psi[j][t] = ind;
        }
    }

    // termination
    TFloat maxval = delta[0][T - 1];
    unsigned int ind = 0;
    for (int i = 1; i < K; ++i) {
        if (delta[i][T - 1] > maxval) {
            maxval = delta[i][T - 1];
            ind = i;
        }
    }

    // traceback
    path[T - 1] = ind;
    for (int t = static_cast<int>(T) - 2; t >= 0; --t) {
        path[t] = psi[path[t + 1]][t + 1];
    }  
}

template<class TFloat>
void SLMSeg<TFloat>::bioviterbii_SLM(
    const std::vector<TFloat>& etav, 
    const std::vector<std::vector<TFloat>>& P, 
    const std::vector<std::vector<TFloat>>& Emission,
    int T, 
    int K, 
    std::vector<unsigned int>& path, 
    std::vector<std::vector<unsigned int>>& psi) 
{
    std::vector<std::vector<TFloat>> delta(K, std::vector<TFloat>(T, 0.0));

    // Initialization
    for (int i = 0; i < K; ++i) {
        delta[i][0] = etav[i] + Emission[i][0];
        psi[i][0] = 0;
    }

    // Recursion
    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < K; ++j) {
            TFloat max_score = delta[0][t - 1] + P[0][j];
            unsigned int max_index = 0;

            for (int i = 1; i < K; ++i) {
                TFloat score = delta[i][t - 1] + P[i][j];
                if (score > max_score) {
                    max_score = score;
                    max_index = i;
                }
            }

            delta[j][t] = max_score + Emission[j][t];
            psi[j][t] = max_index;
        }
    }

    // Termination
    TFloat maxval = delta[0][T - 1];
    unsigned int maxind = 0;
    for (int i = 1; i < K; ++i) {
        if (delta[i][T - 1] > maxval) {
            maxval = delta[i][T - 1];
            maxind = i;
        }
    }

    // Traceback
    path[T - 1] = maxind;
    for (int t = T - 2; t >= 0; --t) {
        path[t] = psi[path[t + 1]][t + 1];
    }
}

template<class TFloat>
void SLMSeg<TFloat>::param_est_seq() 
{
    // smu_: Between-segment variation
    // sepsilon_: Within-segment noise
    data_matrix.emplace_back(signal_data_); // log2R
    int Nexp = data_matrix.size();

    std::vector<TFloat> sigmax(Nexp);

    for (int i = 0; i < Nexp; ++i) {
        // for single sample
        std::vector<TFloat>& row = data_matrix[i]; 
        std::vector<TFloat> centered(row.size());
        // caculate median
        std::vector<TFloat> tmp = row;
        std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end()); // sort nth element
        TFloat median = tmp[tmp.size() / 2]; 
        // caculate absolute deviation from the median
        for (int j = 0; j < row.size(); ++j) {
            centered[j] = std::abs(row[j] - median);
        }
        // caculate mad (median absolute deviation)
        std::nth_element(centered.begin(), centered.begin() + centered.size() / 2, centered.end());
        TFloat mad = centered[centered.size() / 2];
        mad *= 1.4826022185; // consistency factor ???
        sigmax[i] = mad * mad;
        smu_.push_back(std::sqrt(omega_ * sigmax[0]));
        sepsilon_.push_back(std::sqrt((1.0 - omega_) * sigmax[0]));
        mi_.push_back(0.0);
    }
    //std::cout << std::fixed << std::setprecision(6);
    //std::cout << "Estimated parameters:\n";
    //std::cout << "  mi       = " << mi_[0] << "\n";
    //std::cout << "  smu      = " << smu_[0] << "\n";
    //std::cout << "  sepsilon = " << sepsilon_[0] << "\n";
}

template<class TFloat>
void SLMSeg<TFloat>::muk_est() 
{
    int Nexp = data_matrix.size();

    if (Nexp == 1) {
        std::vector<TFloat> temp;
        temp.reserve(21);
        for (TFloat val = -1.0; val <= 1.0 + 1e-8; val += 0.1)
            temp.push_back(val);
        muk_.push_back(temp);
    } // else > 1 ... (to do)

    //std::cout << std::fixed << std::setprecision(6);
    //std::cout << "Estimated Muk:\n";
    //for (size_t i = 0; i < muk_[0].size(); ++i) {
    //    std::cout << "  muk[" << i << "] = " << muk_[0][i] << "\n";
    //}
}

template<class TFloat>
void SLMSeg<TFloat>::joint_seg()
{
    // 檢查資料
    int NExp = data_matrix.size();
    if (NExp == 0 || data_matrix[0].empty())
        throw std::invalid_argument("Data matrix is empty.");

    int T = data_matrix[0].size();           // 序列長度
    int K0 = muk_[0].size();                 // 狀態數（μ的組合數）
    std::vector<TFloat> etav(K0, std::log(1.0 / K0));  // etav = log(1/K0)

    // Step 1: 建立轉移矩陣與機率表
    std::vector<TFloat> G(K0, 0.0);                            // GVECT
    std::vector<std::vector<TFloat>> P(K0, std::vector<TFloat>(K0, 0.0)); // transition matrix
    std::vector<std::vector<TFloat>> Emission(K0, std::vector<TFloat>(T, 0.0)); // emission log likelihood

    // Step 2: 執行 transemisi（同質版本）
    transemisi_SLM(eta_, K0, NExp, T, G, P, Emission);
    std::cerr << "-Transemisi (homogeneous) Completed\n";

    // Step 3: Viterbi 路徑與 psi 表
    std::vector<unsigned int> path(T, 0);
    std::vector<std::vector<unsigned int>> psi(K0, std::vector<unsigned int>(T, 0));
    bioviterbii_SLM(etav, P, Emission, T, K0, path, psi);
    std::cerr << "-bioviterbi Completed\n";

    // Step 4: 取得斷點位置
    total_pred_break_.clear();
    total_pred_break_ = get_breaks(path);
}

template<class TFloat>
void SLMSeg<TFloat>::joint_seg_in()
{
    // CovPos = diff(Pos)
    // CovPosNorm = CovPos / stepeta(dnorm)
    // Pr = etavec = eta + (1 - eta) * exp(log(eta) / CovPosNorm) 
    std::vector<TFloat> etavec(pos_data_.size() - 1);
    for (int i = 0; i < etavec.size(); ++i) {
       TFloat cov_pos = static_cast<TFloat>(pos_data_[i + 1]) - static_cast<TFloat>(pos_data_[i]);
       TFloat cov_pos_norm = cov_pos / stepeta_;
       etavec[i] = eta_ + (1.0 - eta_) * std::exp(std::log(eta_) / cov_pos_norm);
    }
    // 樣本數 NExp = sample size
    int NExp = data_matrix.size();
    if (NExp == 0 || data_matrix[0].empty())
    throw std::invalid_argument("Data matrix is empty.");

    // bins pair NCov = length(etavec)
    int NCov = etavec.size();
    // 狀態數（候選的 mean 組合 μ_k）K0 = ncol(muk)
    int K0 = muk_[0].size();
    // bins number
    int T = data_matrix[0].size();
    // etav = log(rep(1, K0) * (1/K0))
    std::vector<TFloat> etav(K0, std::log(1.0 / K0));

    //initialize matrix
    std::vector<TFloat> G(K0, 0.0);
    std::vector<std::vector<TFloat>> P(K0, std::vector<TFloat>(K0 * NCov, 0.0));
    std::vector<std::vector<TFloat>> Emission(K0, std::vector<TFloat>(T, 0.0));
    // do 
    transemisi_HSLM(etavec, NCov, K0, NExp, T, G, P, Emission);
    std::cerr << "-Transemisi Completed \n";

    // initialize matrix
    std::vector<unsigned int> path(T, 0);
    std::vector<std::vector<unsigned int>> psi(K0, std::vector<unsigned int>(T, 0));
    // do 
    bioviterbii_HSLM(etav, P, Emission, T, K0, path, psi);
    std::cerr << "-bioviterbii Completed \n";

    total_pred_break_.clear();
    total_pred_break_ = get_breaks(path);
    //write_breakpoints_to_file("file2.txt");
    //for(int i = 0; i< total_pred_break_.size(); ++i) {
    //    std::cout << total_pred_break_[i] << std::endl;
    //}
}

template<class TFloat>
void SLMSeg<TFloat>::filter_seg() {
    total_pred_break_filtered_.clear(); 
    int n = total_pred_break_.size();
    if (n == 0) return;

    std::vector<unsigned int> indF;
    for (size_t i = 0; i < n - 1; ++i) {
        TFloat controllength = total_pred_break_[i + 1] - total_pred_break_[i];
        if (controllength <= fw_) {
            indF.push_back(i);
        }
    }

    std::set<unsigned int> skip_indices;
    if (!indF.empty()) {
        if (indF[0] == 1)
            indF[0] = 2;

        indF.erase(std::unique(indF.begin(), indF.end()), indF.end());
        skip_indices.insert(indF.begin(), indF.end());
    }

    for (size_t i = 0; i < n; ++i) {
        if (!skip_indices.count(i)) {
            total_pred_break_filtered_.push_back(total_pred_break_[i]);
        }
    }
}    

template<class TFloat>
void SLMSeg<TFloat>::seg_results()
{
    int NExp = data_matrix.size();
    int T = data_matrix[0].size();

    for (int j = 0; j < NExp; ++j) {
        std::vector<TFloat> s(T, 0.0);
        for (int i = 0; i + 1 < total_pred_break_filtered_.size(); ++i) {
            unsigned int start = total_pred_break_filtered_[i];
            unsigned int end = total_pred_break_filtered_[i + 1];

            std::vector<TFloat> segment(data_matrix[j].begin() + start, data_matrix[j].begin() + end);
            std::nth_element(segment.begin(), segment.begin() + segment.size() / 2, segment.end());
            TFloat median = segment[segment.size() / 2];

            for (int t = start; t < end; ++t)
                s[t] = median;
        }
        data_seg_.push_back(std::move(s));
    }
}

template<class TFloat>
void SLMSeg<TFloat>::HSLM() 
{
    std::cout << "Parameter estimation...\n";
    param_est_seq();

    std::cout << "Muk estimation...\n";
    muk_est();

    std::cout << "Joint segmentation...\n";
    joint_seg_in();

    std::cout << "Filtering predictions...\n";
    filter_seg();
    //write_breakpoints_to_file("file2.txt");

    std::cout << "Producing results...\n";
    seg_results();
}

template<class TFloat>
void SLMSeg<TFloat>::SLM() 
{
    std::cout << "Parameter estimation...\n";
    param_est_seq();

    std::cout << "Muk estimation...\n";
    muk_est();

    std::cout << "Joint segmentation...\n";
    joint_seg();

    std::cout << "Filtering predictions...\n";
    filter_seg();
    //write_breakpoints_to_file("file2.txt");

    std::cout << "Producing results...\n";
    seg_results();
}

#endif // SLMSEG_H