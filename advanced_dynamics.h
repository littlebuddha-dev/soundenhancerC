// ./advanced_dynamics.h
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>

// マルチバンドコンプレッサー
class MultibandCompressor {
public:
    struct Band {
        double freq_low, freq_high;
        double threshold_db, ratio, attack_ms, release_ms;
        double makeup_gain_db;
        bool enabled;
        
        // 内部状態
        double envelope = 0.0;
        double attack_coeff = 0.0;
        double release_coeff = 0.0;
        SimpleBiquad lpf, hpf, bpf;
    };
    
    void setup(double sr, const std::vector<Band>& bands) {
        sample_rate_ = sr;
        bands_ = bands;
        
        for (auto& band : bands_) {
            if (!band.enabled) continue;
            
            // フィルター設定
            band.lpf.set_lpf(sr, band.freq_high, 0.707);
            band.hpf.set_hpf(sr, band.freq_low, 0.707);
            
            // エンベロープ係数計算
            band.attack_coeff = std::exp(-1.0 / (band.attack_ms / 1000.0 * sr));
            band.release_coeff = std::exp(-1.0 / (band.release_ms / 1000.0 * sr));
        }
        
        // クロスオーバーネットワーク初期化
        setupCrossoverNetwork();
    }
    
    float process(float input) {
        std::vector<float> band_signals = splitToBands(input);
        std::vector<float> compressed_bands;
        
        for (size_t i = 0; i < bands_.size(); ++i) {
            if (!bands_[i].enabled) {
                compressed_bands.push_back(band_signals[i]);
                continue;
            }
            
            float compressed = compressBand(band_signals[i], bands_[i]);
            compressed_bands.push_back(compressed);
        }
        
        return sumBands(compressed_bands);
    }

private:
    double sample_rate_ = 44100.0;
    std::vector<Band> bands_;
    std::vector<SimpleBiquad> crossover_filters_;
    
    void setupCrossoverNetwork() {
        // リンクウィッツ・ライリー4次クロスオーバー実装
        crossover_filters_.clear();
        for (size_t i = 0; i < bands_.size() - 1; ++i) {
            SimpleBiquad lpf1, lpf2, hpf1, hpf2;
            double freq = bands_[i].freq_high;
            
            lpf1.set_lpf(sample_rate_, freq, 0.707);
            lpf2.set_lpf(sample_rate_, freq, 0.707);
            hpf1.set_hpf(sample_rate_, freq, 0.707);
            hpf2.set_hpf(sample_rate_, freq, 0.707);
            
            crossover_filters_.push_back(lpf1);
            crossover_filters_.push_back(lpf2);
            crossover_filters_.push_back(hpf1);
            crossover_filters_.push_back(hpf2);
        }
    }
    
    std::vector<float> splitToBands(float input) {
        std::vector<float> bands(bands_.size(), 0.0f);
        
        if (bands_.size() == 1) {
            bands[0] = input;
            return bands;
        }
        
        // 簡略化された帯域分割（実際にはより複雑なクロスオーバーが必要）
        float low = input;
        for (size_t i = 0; i < bands_.size(); ++i) {
            if (i == 0) {
                bands[i] = bands_[i].lpf.process(low);
            } else if (i == bands_.size() - 1) {
                bands[i] = bands_[i].hpf.process(low);
            } else {
                float temp = bands_[i].hpf.process(low);
                bands[i] = bands_[i].lpf.process(temp);
            }
        }
        
        return bands;
    }
    
    float compressBand(float input, Band& band) {
        float abs_input = std::abs(input);
        
        // エンベロープフォロワー
        if (abs_input > band.envelope) {
            band.envelope = band.attack_coeff * band.envelope + (1.0 - band.attack_coeff) * abs_input;
        } else {
            band.envelope = band.release_coeff * band.envelope + (1.0 - band.release_coeff) * abs_input;
        }
        
        // 圧縮計算
        double threshold_linear = db_to_linear(band.threshold_db);
        double gain_reduction = 1.0;
        
        if (band.envelope > threshold_linear) {
            double over_threshold_db = 20.0 * std::log10(band.envelope / threshold_linear);
            double compressed_db = over_threshold_db / band.ratio;
            gain_reduction = db_to_linear(compressed_db - over_threshold_db);
        }
        
        // メイクアップゲイン適用
        double makeup_gain = db_to_linear(band.makeup_gain_db);
        
        return input * gain_reduction * makeup_gain;
    }
    
    float sumBands(const std::vector<float>& bands) {
        float sum = 0.0f;
        for (float band : bands) {
            sum += band;
        }
        return sum;
    }
};

// アナログ風サチュレーション
class AnalogSaturation {
public:
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        drive_ = params.value("drive", 1.0);
        mix_ = params.value("mix", 0.3);
        type_ = params.value("type", "tube"); // "tube", "tape", "transformer"
        
        // DCブロッキングフィルター
        dc_blocker_.set_hpf(sr, 5.0, 0.707);
        
        // アンチエイリアシングフィルター
        anti_alias_.set_lpf(sr, sr * 0.4, 0.707);
    }
    
    float process(float input) {
        float processed = input;
        
        // プリフィルタリング
        processed = anti_alias_.process(processed);
        
        // サチュレーション
        if (type_ == "tube") {
            processed = tubeSaturation(processed);
        } else if (type_ == "tape") {
            processed = tapeSaturation(processed);
        } else if (type_ == "transformer") {
            processed = transformerSaturation(processed);
        }
        
        // ポストフィルタリング
        processed = dc_blocker_.process(processed);
        
        // ミックス
        return input * (1.0f - mix_) + processed * mix_;
    }

private:
    double sample_rate_ = 44100.0;
    double drive_ = 1.0;
    double mix_ = 0.3;
    std::string type_ = "tube";
    SimpleBiquad dc_blocker_, anti_alias_;
    
    float tubeSaturation(float x) {
        x *= drive_;
        // 真空管の非線形特性をモデル化
        float sign = (x >= 0.0f) ? 1.0f : -1.0f;
        float abs_x = std::abs(x);
        
        if (abs_x < 0.3f) {
            return x * (1.0f + 0.1f * abs_x);
        } else if (abs_x < 0.7f) {
            return sign * (0.3f + 0.4f * std::tanh((abs_x - 0.3f) * 2.5f));
        } else {
            return sign * std::tanh(abs_x * 0.8f);
        }
    }
    
    float tapeSaturation(float x) {
        x *= drive_;
        // テープの非線形特性
        return std::tanh(x * 1.2f) * 0.9f;
    }
    
    float transformerSaturation(float x) {
        x *= drive_;
        // トランスフォーマーの非線形特性
        float x2 = x * x;
        return x * (1.0f - x2 * 0.15f + x2 * x2 * 0.02f);
    }
};