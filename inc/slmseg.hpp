#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numbers>
#include <set>
#include <stdexcept>
#include <execution>
#include <array>
#include <algorithm> 
#include <unordered_set>

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

    void transemisi_HSLM(const std::vector<TFloat>& etavec, size_t NCov, size_t K, size_t NExp, size_t T,
                    std::vector<TFloat>& G, std::vector<std::vector<TFloat>>& P, std::vector<std::vector<TFloat>>& Emission);

    void transemisi_SLM(const TFloat eta, size_t K, size_t NExp, size_t T,
                    std::vector<TFloat>& G, std::vector<std::vector<TFloat>>& P, std::vector<std::vector<TFloat>>& Emission);                  
                    
    void bioviterbii_HSLM(const std::vector<TFloat>& etav, const std::vector<std::vector<TFloat>>& P, 
                    const std::vector<std::vector<TFloat>>& Emission, size_t  T, size_t  K, std::vector<unsigned int>& path, 
                    std::vector<std::vector<unsigned int>>& psi);
    
    void bioviterbii_SLM(const std::vector<TFloat>& etav, const std::vector<std::vector<TFloat>>& P, 
                    const std::vector<std::vector<TFloat>>& Emission, size_t  T, size_t  K, std::vector<unsigned int>& path, 
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

    static inline constexpr std::array<TFloat, 11> BIN_VEC = {
        TFloat(-1.0), TFloat(-0.8), TFloat(-0.6), TFloat(-0.4), TFloat(-0.2),
        TFloat(0.0),
        TFloat(0.2), TFloat(0.4), TFloat(0.6), TFloat(0.8), TFloat(1.0)
    };

    static inline constexpr std::array<TFloat, 10> BIN_MEAN = {
        TFloat(-0.9), TFloat(-0.7), TFloat(-0.5), TFloat(-0.3),
        TFloat(0.0), TFloat(0.0),
        TFloat(0.3), TFloat(0.5), TFloat(0.7), TFloat(0.9)
    };  
    
    static inline constexpr std::array<TFloat, 21> MUK_SINGLE = {
        TFloat(-1.0), TFloat(-0.9), TFloat(-0.8), TFloat(-0.7), TFloat(-0.6),
        TFloat(-0.5), TFloat(-0.4), TFloat(-0.3), TFloat(-0.2), TFloat(-0.1),
        TFloat(0.0),  TFloat(0.1),  TFloat(0.2),  TFloat(0.3),  TFloat(0.4),
        TFloat(0.5),  TFloat(0.6),  TFloat(0.7),  TFloat(0.8),  TFloat(0.9),
        TFloat(1.0)
    };

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

    size_t T = path.size();
    unsigned int last_state = path[0];
    breakpoints.push_back(0);  // 起點一定是斷點

    for (size_t t = 1; t < T; ++t) {
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
    size_t NCov,
    size_t K,
    size_t NExp,
    size_t T,
    std::vector<TFloat>& G, 
    std::vector<std::vector<TFloat>>& P, 
    std::vector<std::vector<TFloat>>& Emission) 
{
    TFloat PI = std::numbers::pi;
    // Step 1: GVECT (a simplify gaussion log-likelihood)
    std::vector<TFloat> gvect(K, 0.0);
    for (size_t i = 0; i < K; ++i) {
        TFloat gsum = 0.0;
        for (size_t j = 0; j < NExp; ++j) {
            TFloat diff = muk_[j][i] - mi_[j];
            gsum += - (diff * diff) / (2 * smu_[j] * smu_[j]);
        }
        gvect[i] = gsum;
    }

    // Step 2: normalize GVECT (likelihood) to G (probability)
    TFloat norm = gvect[0];
    for (size_t i = 1; i < K; ++i)
        norm = elnsum(norm, gvect[i]);
    for (size_t i = 0; i < K; ++i)
        G[i] = gvect[i] - norm;

    // Step 3: transition matrix P [NCov][K][K] → flattened to [K * NCov][K]
    for (size_t cov = 0; cov < NCov; ++cov) {
        for (size_t i = 0; i < K; ++i) {
            for (size_t j = 0; j < K; ++j) {
                size_t row = cov * K + i;
                if (i == j)
                    P[row][j] = elnsum(std::log(1 - etavec[cov]), std::log(etavec[cov]) + G[j]);
                else
                    P[row][j] = std::log(etavec[cov]) + G[j];
            }
        }
    }

    // Step 4: emission log-likelihood [T][K] 
    for (size_t t = 0; t < T; ++t) {
        for (size_t k = 0; k < K; ++k) {
            TFloat val = 0.0;
            for (size_t j = 0; j < NExp; ++j) {
                TFloat diff = data_matrix[j][t] - muk_[j][k];
                val += -0.5 * (diff * diff) / (sepsilon_[j] * sepsilon_[j]) -
                    std::log(std::sqrt(2.0 * PI) * sepsilon_[j]);
            }
            Emission[t][k] = val;
        }
    }

}

template<class TFloat>
void SLMSeg<TFloat>::transemisi_SLM(
    const TFloat eta,
    size_t K,
    size_t NExp,
    size_t T,
    std::vector<TFloat>& G, 
    std::vector<std::vector<TFloat>>& P, 
    std::vector<std::vector<TFloat>>& Emission) 
{
    TFloat PI = std::numbers::pi;

    // Step 1: GVECT
    std::vector<TFloat> gvect(K, 0.0);
    for (size_t i = 0; i < K; ++i) {
        TFloat gsum = 0.0;
        for (size_t j = 0; j < NExp; ++j) {
            TFloat diff = muk_[j][i] - mi_[j];
            gsum += - (diff * diff) / (2 * smu_[j] * smu_[j]);
        }
        gvect[i] = gsum;
    }

    // Step 2: Normalize GVECT to G
    TFloat norm = gvect[0];
    for (size_t i = 1; i < K; ++i)
        norm = elnsum(norm, gvect[i]);
    for (size_t i = 0; i < K; ++i)
        G[i] = gvect[i] - norm;

    // Step 3: Transition matrix P (K x K)
    for (size_t i = 0; i < K; ++i) {
        for (size_t j = 0; j < K; ++j) {
            if (i == j)
                P[i][j] = elnsum(std::log(1.0 - eta), std::log(eta) + G[j]);
            else
                P[i][j] = std::log(eta) + G[j];
        }
    }

    // Step 4: emission log-likelihood [T][K] 
    for (size_t t = 0; t < T; ++t) {
        for (size_t k = 0; k < K; ++k) {
            TFloat val = 0.0;
            for (size_t j = 0; j < NExp; ++j) {
                TFloat diff = data_matrix[j][t] - muk_[j][k];
                val += -0.5 * (diff * diff) / (sepsilon_[j] * sepsilon_[j]) -
                    std::log(std::sqrt(2.0 * PI) * sepsilon_[j]);
            }
            Emission[t][k] = val;
        }
    }
}

template<class TFloat>
void SLMSeg<TFloat>::bioviterbii_HSLM(
    const std::vector<TFloat>& etav, 
    const std::vector<std::vector<TFloat>>& P, 
    const std::vector<std::vector<TFloat>>& Emission,
    size_t T, 
    size_t K, 
    std::vector<unsigned int>& path, 
    std::vector<std::vector<unsigned int>>& psi) 
{
    std::vector<std::vector<TFloat>> delta(T, std::vector<TFloat>(K, 0.0));
    
    // initialization
    for (size_t i = 0; i < K; ++i) {
        delta[0][i] = etav[i] + Emission[0][i];
        psi[0][i] = 0;
    }

    // recursion
    for (size_t t = 1; t < T; ++t) {
        for (size_t j = 0; j < K; ++j) {
            size_t offset = (t - 1) * K; // offset + i, i=0 
            TFloat nummax = delta[t - 1][0] + P[offset][j]; // P[0][j][t]
            unsigned int ind = 0;

            for (size_t i = 1; i < K; ++i) {
                TFloat score = delta[t - 1][i] + P[offset + i][j];// P[i][j][t-1]
                if (score > nummax) {
                    nummax = score;
                    ind = i;
                }
            }

            delta[t][j] = nummax + Emission[t][j];
            psi[t][j] = ind;
        }
    }

    // termination
    TFloat maxval = delta[T - 1][0];
    unsigned int ind = 0;
    for (size_t i = 1; i < K; ++i) {
        if (delta[T - 1][i] > maxval) {
            maxval = delta[T - 1][i];
            ind = i;
        }
    }

    // traceback
    path[T - 1] = ind;
    for (int t = static_cast<int>(T) - 2; t >= 0; --t) {
        path[t] = psi[t + 1][path[t + 1]];
    }  
}

template<class TFloat>
void SLMSeg<TFloat>::bioviterbii_SLM(
    const std::vector<TFloat>& etav, 
    const std::vector<std::vector<TFloat>>& P, 
    const std::vector<std::vector<TFloat>>& Emission,
    size_t T, 
    size_t K, 
    std::vector<unsigned int>& path, 
    std::vector<std::vector<unsigned int>>& psi) 
{
    std::vector<std::vector<TFloat>> delta(T, std::vector<TFloat>(K, 0.0));
    
    // initialization
    for (size_t i = 0; i < K; ++i) {
        delta[0][i] = etav[i] + Emission[0][i];
        psi[0][i] = 0;
    }

    // Recursion
    for (size_t t = 1; t < T; ++t) {
        for (size_t j = 0; j < K; ++j) {
            TFloat nummax = delta[t - 1][0] + P[0][j];
            unsigned int ind = 0;

            for (size_t i = 1; i < K; ++i) {
                TFloat score = delta[t - 1][i] + P[i][j];
                if (score > nummax) {
                    nummax = score;
                    ind = i;
                }
            }

            delta[t][j] = nummax + Emission[t][j];
            psi[t][j] = ind;
        }
    }

    // Termination
    TFloat maxval = delta[T - 1][0];
    unsigned int ind = 0;
    for (size_t i = 1; i < K; ++i) {
        if (delta[T - 1][i] > maxval) {
            maxval = delta[T - 1][i];
            ind = i;
        }
    }

    // traceback
    path[T - 1] = ind;
    for (int t = static_cast<int>(T) - 2; t >= 0; --t) {
        path[t] = psi[t + 1][path[t + 1]];
    }  
}

template<class TFloat>
void SLMSeg<TFloat>::param_est_seq() 
{   
    data_matrix.clear();
    smu_.clear();
    sepsilon_.clear();
    mi_.clear();

    // smu_: Between-segment variation
    // sepsilon_: Within-segment noise
    data_matrix.emplace_back(signal_data_); // log2R
    size_t Nexp = data_matrix.size();
    if (Nexp < 1) {
        std::cerr << "Must at least one sample" <<std::endl;
        std::exit(1); 
    }

    std::vector<TFloat> sigmax(Nexp);

    for (size_t i = 0; i < Nexp; ++i) {
        // for single sample
        std::vector<TFloat>& row = data_matrix[i]; 
        std::vector<TFloat> centered(row.size());
        // caculate median
        std::vector<TFloat> tmp = row;
        std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end()); // sort nth element
        TFloat median = tmp[tmp.size() / 2]; 
        // caculate absolute deviation from the median
        for (size_t j = 0; j < row.size(); ++j) {
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
    muk_.clear();
    size_t Nexp = data_matrix.size();

    if (Nexp == 1) {
        muk_.emplace_back(MUK_SINGLE.begin(), MUK_SINGLE.end());
        return;
    }

    // Clip to [-1, 1]
    for (auto& sample : data_matrix) {
        for (auto& val : sample) {
            val = std::max(TFloat(-1.0), std::min(val, TFloat(1.0)));
        }
    }

    // 1. Moving average smoothing
    std::vector<std::vector<TFloat>> data_matrix_smooth(Nexp);

    for (size_t i = 0; i < Nexp; ++i) {
        const auto& sample = data_matrix[i];
        size_t len = sample.size();
        std::vector<TFloat> smoothed;

        for (size_t j = 0; j + mw_ <= len; ++j) {
            TFloat sum = 0;
            for (size_t k = 0; k < mw_; ++k) {
                sum += sample[j + k];
            }
            smoothed.push_back(sum / mw_);
        }

        data_matrix_smooth[i] = std::move(smoothed);
    }

    // 2. Quantize with BIN_VEC + BIN_MEAN
    for (size_t i = 0; i < Nexp; ++i) {
        for (size_t t = 0; t < data_matrix_smooth[i].size(); ++t) {
            TFloat& val = data_matrix_smooth[i][t];
            auto it = std::upper_bound(BIN_VEC.begin(), BIN_VEC.end(), val);
            if (it == BIN_VEC.begin()) {
                val = BIN_MEAN.front();
            } else if (it == BIN_VEC.end()) {
                val = BIN_MEAN.back();
            } else {
                size_t idx = std::distance(BIN_VEC.begin(), it) - 1;
                val = BIN_MEAN[idx];
            }
        }
    }

    // 3. Unique column extraction, skip first column (col=1)
    size_t Tprime = data_matrix_smooth[0].size();
    std::unordered_set<std::string> seen;
    for (size_t col = 0; col < Tprime; ++col) {
        std::vector<TFloat> col_vector(Nexp);
        for (size_t row = 0; row < Nexp; ++row)
            col_vector[row] = data_matrix_smooth[row][col];

        std::ostringstream oss;
        for (TFloat v : col_vector) oss << v << ',';
        std::string key = oss.str();

        if (seen.insert(key).second) {
            muk_.emplace_back(std::move(col_vector));
        }
    }
}

template<class TFloat>
void SLMSeg<TFloat>::joint_seg()
{
    // 檢查資料
    size_t NExp = data_matrix.size();
    if (NExp == 0 || data_matrix[0].empty())
        throw std::invalid_argument("Data matrix is empty.");

    size_t T = data_matrix[0].size();           // 序列長度
    size_t K0 = muk_[0].size();                 // 狀態數（μ的組合數）
    std::vector<TFloat> etav(K0, std::log(1.0 / K0));  // etav = log(1/K0)

    // Step 1: 建立轉移矩陣與機率表
    std::vector<TFloat> G(K0, 0.0);                            // GVECT
    std::vector<std::vector<TFloat>> P(K0, std::vector<TFloat>(K0, 0.0)); // transition matrix
    std::vector<std::vector<TFloat>> Emission(T, std::vector<TFloat>(K0, 0.0)); // emission log likelihood

    // Step 2: 執行 transemisi（同質版本）
    transemisi_SLM(eta_, K0, NExp, T, G, P, Emission);
    std::cerr << "-Transemisi (homogeneous) Completed\n";

    // Step 3: Viterbi 路徑與 psi 表
    std::vector<unsigned int> path(T, 0);
    std::vector<std::vector<unsigned int>> psi(T, std::vector<unsigned int>(K0, 0));
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
    for (size_t i = 0; i < etavec.size(); ++i) {
       TFloat cov_pos = static_cast<TFloat>(pos_data_[i + 1]) - static_cast<TFloat>(pos_data_[i]);
       TFloat cov_pos_norm = cov_pos / stepeta_;
       etavec[i] = eta_ + (1.0 - eta_) * std::exp(std::log(eta_) / cov_pos_norm);
    }
    // 樣本數 NExp = sample size
    size_t NExp = data_matrix.size();
    if (NExp == 0 || data_matrix[0].empty())
    throw std::invalid_argument("Data matrix is empty.");

    // bins pair NCov = length(etavec)
    size_t NCov = etavec.size();
    // 狀態數（候選的 mean 組合 μ_k）K0 = ncol(muk)
    size_t K0 = muk_[0].size();
    // bins number
    size_t T = data_matrix[0].size();
    // etav = log(rep(1, K0) * (1/K0))
    std::vector<TFloat> etav(K0, std::log(1.0 / K0));

    //initialize matrix
    std::vector<TFloat> G(K0, 0.0);
    std::vector<std::vector<TFloat>> P(K0 * NCov, std::vector<TFloat>(K0, 0.0));
    std::vector<std::vector<TFloat>> Emission(T, std::vector<TFloat>(K0, 0.0));
    // do 
    transemisi_HSLM(etavec, NCov, K0, NExp, T, G, P, Emission);
    std::cerr << "-Transemisi Completed \n";

    // initialize matrix
    std::vector<unsigned int> path(T, 0);
    std::vector<std::vector<unsigned int>> psi(T, std::vector<unsigned int>(K0, 0));
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
    size_t n = total_pred_break_.size();
    if (n <= 0) return;

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
    size_t NExp = data_matrix.size();
    size_t T = data_matrix[0].size();

    for (size_t j = 0; j < NExp; ++j) {
        std::vector<TFloat> s(T, 0.0);
        for (size_t i = 0; i + 1 < total_pred_break_filtered_.size(); ++i) {
            unsigned int start = total_pred_break_filtered_[i];
            unsigned int end = total_pred_break_filtered_[i + 1];

            std::vector<TFloat> segment(data_matrix[j].begin() + start, data_matrix[j].begin() + end);
            std::nth_element(segment.begin(), segment.begin() + segment.size() / 2, segment.end());
            TFloat median = segment[segment.size() / 2];

            for (size_t t = start; t < end; ++t)
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
