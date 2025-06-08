// ./advanced_eq_harmonics.h
#pragma once
#include <vector>
#include <cmath>
#include <complex>

// 線形位相EQ（FFTベース）
class LinearPhaseEQ {
public:
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        fft_size_ = params.value("fft_size", 2048);
        overlap_ = params.value("overlap", 0.75);
        enabled_ = params.value("enabled", true);
        
        // EQカーブ設定
        setupEQCurve(params.value("eq_bands", json::array()));
        
        // オーバーラップアド処理用バッファ
        hop_size_ = static_cast<size_t>(fft_size_ * (1.0 - overlap_));
        input_buffer_.resize(fft_size_, 0.0f);
        output_buffer_.resize(fft_size_, 0.0f);
        overlap_buffer_.resize(fft_size_, 0.0f);
        
        // ウィンドウ関数（ハニング窓）
        window_.resize(fft_size_);
        for (size_t i = 0; i < fft_size_; ++i) {
            window_[i] = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (fft_size_ - 1)));
        }
    }
    
    std::vector<float> process(const std::vector<float>& input) {
        if (!enabled_) return input;
        
        std::vector<float> output;
        output.reserve(input.size());
        
        for (float sample : input) {
            // バッファに蓄積
            input_buffer_[buffer_pos_] = sample;
            buffer_pos_ = (buffer_pos_ + 1) % fft_size_;
            
            samples_since_last_process_++;
            
            // hop_size_ごとに処理
            if (samples_since_last_process_ >= hop_size_) {
                processBlock();
                samples_since_last_process_ = 0;
            }
            
            // 出力
            output.push_back(output_buffer_[output_pos_]);
            output_pos_ = (output_pos_ + 1) % fft_size_;
        }
        
        return output;
    }

private:
    double sample_rate_ = 44100.0;
    size_t fft_size_ = 2048;
    double overlap_ = 0.75;
    size_t hop_size_;
    bool enabled_ = true;
    
    std::vector<float> input_buffer_, output_buffer_, overlap_buffer_;
    std::vector<float> window_;
    std::vector<std::complex<float>> eq_response_;
    
    size_t buffer_pos_ = 0;
    size_t output_pos_ = 0;
    size_t samples_since_last_process_ = 0;
    
    void setupEQCurve(const json& bands) {
        eq_response_.resize(fft_size_ / 2 + 1);
        
        // 初期化（フラット）
        for (auto& response : eq_response_) {
            response = std::complex<float>(1.0f, 0.0f);
        }
        
        // 各バンドを処理
        for (const auto& band : bands) {
            if (!band.value("enabled", true)) continue;
            
            double freq = band.value("freq", 1000.0);
            double gain_db = band.value("gain_db", 0.0);
            double q = band.value("q", 1.0);
            std::string type = band.value("type", "peaking");
            
            applyEQBand(freq, gain_db, q, type);
        }
    }
    
    void applyEQBand(double freq, double gain_db, double q, const std::string& type) {
        double gain_linear = db_to_linear(gain_db);
        
        for (size_t bin = 0; bin < eq_response_.size(); ++bin) {
            double bin_freq = bin * sample_rate_ / (2.0 * (eq_response_.size() - 1));
            double response = 1.0;
            
            if (type == "peaking") {
                double ratio = bin_freq / freq;
                double bandwidth = freq / q;
                response = 1.0 + (gain_linear - 1.0) / 
                          (1.0 + std::pow((bin_freq - freq) / (bandwidth / 2.0), 2.0));
            } else if (type == "lowpass") {
                double ratio = bin_freq / freq;
                response = 1.0 / std::sqrt(1.0 + std::pow(ratio * q, 2.0));
                if (gain_db != 0.0) response *= gain_linear;
            } else if (type == "highpass") {
                double ratio = freq / bin_freq;
                response = 1.0 / std::sqrt(1.0 + std::pow(ratio * q, 2.0));
                if (gain_db != 0.0) response *= gain_linear;
            }
            
            eq_response_[bin] *= std::complex<float>(response, 0.0f);
        }
    }
    
    void processBlock() {
        // 実装は簡略化（実際にはFFT/IFFT処理が必要）
        // ここでは概念的な処理のみ示す
    }
};

// ハーモニックエンハンサー
class HarmonicEnhancer {
public:
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        enabled_ = params.value("enabled", true);
        drive_ = params.value("drive", 0.3);
        even_harmonics_ = params.value("even_harmonics", 0.2);
        odd_harmonics_ = params.value("odd_harmonics", 0.3);
        mix_ = params.value("mix", 0.25);
        
        // バンドパスフィルター（処理対象周波数帯域）
        bandpass_.set_hpf(sr, params.value("freq_low", 100.0), 0.707);
        lowpass_.set_lpf(sr, params.value("freq_high", 8000.0), 0.707);
        
        // DCブロッキング
        dc_blocker_.set_hpf(sr, 10.0, 0.707);
    }
    
    float process(float input) {
        if (!enabled_) return input;
        
        // 処理対象帯域を抽出
        float filtered = bandpass_.process(input);
        filtered = lowpass_.process(filtered);
        
        // ハーモニック生成
        float harmonics = generateHarmonics(filtered);
        
        // DCブロッキング
        harmonics = dc_blocker_.process(harmonics);
        
        // ミックス
        return input * (1.0f - mix_) + (input + harmonics) * mix_;
    }

private:
    double sample_rate_ = 44100.0;
    bool enabled_ = true;
    double drive_ = 0.3;
    double even_harmonics_ = 0.2;
    double odd_harmonics_ = 0.3;
    double mix_ = 0.25;
    
    SimpleBiquad bandpass_, lowpass_, dc_blocker_;
    
    float generateHarmonics(float input) {
        float enhanced = input * drive_;
        
        // 偶次高調波（ウォーム感）
        float even = enhanced * enhanced * even_harmonics_;
        
        // 奇次高調波（明瞭感）
        float odd = std::tanh(enhanced * 2.0f) * odd_harmonics_;
        
        // チューブ風の非線形性
        float tube_like = enhanced / (1.0f + std::abs(enhanced) * 0.3f) * 0.1f;
        
        return even + odd + tube_like;
    }
};

// スペクトラルゲート（ノイズ除去）
class SpectralGate {
public:
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        enabled_ = params.value("enabled", false);
        threshold_db_ = params.value("threshold_db", -60.0);
        attack_ms_ = params.value("attack_ms", 5.0);
        release_ms_ = params.value("release_ms", 100.0);
        
        // スムージングフィルター
        double attack_freq = 1000.0 / (2.0 * M_PI * attack_ms_);
        double release_freq = 1000.0 / (2.0 * M_PI * release_ms_);
        
        attack_smoother_.set_lpf(sr, attack_freq, 0.707);
        release_smoother_.set_lpf(sr, release_freq, 0.707);
    }
    
    float process(float input) {
        if (!enabled_) return input;
        
        float abs_input = std::abs(input);
        float input_db = (abs_input > 0.0f) ? 20.0f * std::log10f(abs_input) : -120.0f;
        
        // ゲート判定
        bool gate_open = input_db > threshold_db_;
        float target_gain = gate_open ? 1.0f : 0.0f;
        
        // スムージング
        if (target_gain > current_gain_) {
            current_gain_ = attack_smoother_.process(target_gain);
        } else {
            current_gain_ = release_smoother_.process(target_gain);
        }
        
        return input * current_gain_;
    }

private:
    double sample_rate_ = 44100.0;
    bool enabled_ = false;
    double threshold_db_ = -60.0;
    double attack_ms_ = 5.0;
    double release_ms_ = 100.0;
    
    float current_gain_ = 0.0f;
    SimpleBiquad attack_smoother_, release_smoother_;
};