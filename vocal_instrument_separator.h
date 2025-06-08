// ./vocal_instrument_separator.h
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <array>

// M/S (Mid-Side) ベースの分離プロセッサー
class MSVocalInstrumentSeparator {
public:
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        enabled_ = params.value("enabled", true);
        
        // ボーカル強化設定
        vocal_enhance_ = params.value("vocal_enhance", 0.3);
        vocal_center_freq_ = params.value("vocal_center_freq", 2500.0);
        vocal_bandwidth_ = params.value("vocal_bandwidth", 2000.0);
        
        // 楽器強化設定
        instrument_enhance_ = params.value("instrument_enhance", 0.2);
        stereo_width_ = params.value("stereo_width", 1.2);
        
        // バンドパスフィルター（ボーカル帯域）
        vocal_bandpass_low_.set_hpf(sr, vocal_center_freq_ - vocal_bandwidth_/2, 0.707);
        vocal_bandpass_high_.set_lpf(sr, vocal_center_freq_ + vocal_bandwidth_/2, 0.707);
        
        // 楽器用ローパス・ハイパス
        instrument_low_.set_lpf(sr, 800.0, 0.8);  // 低域楽器
        instrument_high_.set_hpf(sr, 6000.0, 0.8); // 高域楽器
        
        // ダイナミック検出用エンベロープフォロワー
        setupEnvelopeFollowers(sr);
    }
    
    std::pair<float, float> process(float left, float right) {
        if (!enabled_) return {left, right};
        
        // M/S変換
        float mid = (left + right) * 0.5f;
        float side = (left - right) * 0.5f;
        
        // ボーカル・楽器成分の検出と分離
        auto separated = detectAndSeparate(mid, side);
        float enhanced_mid = separated.first;
        float enhanced_side = separated.second;
        
        // L/R復元
        float enhanced_left = enhanced_mid + enhanced_side;
        float enhanced_right = enhanced_mid - enhanced_side;
        
        return {enhanced_left, enhanced_right};
    }

private:
    double sample_rate_ = 44100.0;
    bool enabled_ = true;
    double vocal_enhance_ = 0.3;
    double vocal_center_freq_ = 2500.0;
    double vocal_bandwidth_ = 2000.0;
    double instrument_enhance_ = 0.2;
    double stereo_width_ = 1.2;
    
    SimpleBiquad vocal_bandpass_low_, vocal_bandpass_high_;
    SimpleBiquad instrument_low_, instrument_high_;
    
    // エンベロープフォロワー
    float vocal_envelope_ = 0.0f;
    float instrument_envelope_ = 0.0f;
    float vocal_attack_coeff_ = 0.0f;
    float vocal_release_coeff_ = 0.0f;
    float inst_attack_coeff_ = 0.0f;
    float inst_release_coeff_ = 0.0f;
    
    void setupEnvelopeFollowers(double sr) {
        // ボーカル検出用（高速応答）
        vocal_attack_coeff_ = std::exp(-1.0 / (0.003 * sr));   // 3ms
        vocal_release_coeff_ = std::exp(-1.0 / (0.1 * sr));    // 100ms
        
        // 楽器検出用（中速応答）
        inst_attack_coeff_ = std::exp(-1.0 / (0.01 * sr));     // 10ms
        inst_release_coeff_ = std::exp(-1.0 / (0.05 * sr));    // 50ms
    }
    
    std::pair<float, float> detectAndSeparate(float mid, float side) {
        // ボーカル成分検出（主にMidチャンネル）
        float vocal_signal = vocal_bandpass_high_.process(
            vocal_bandpass_low_.process(mid)
        );
        float vocal_level = std::abs(vocal_signal);
        
        // エンベロープフォロワー更新
        if (vocal_level > vocal_envelope_) {
            vocal_envelope_ = vocal_attack_coeff_ * vocal_envelope_ + 
                             (1.0f - vocal_attack_coeff_) * vocal_level;
        } else {
            vocal_envelope_ = vocal_release_coeff_ * vocal_envelope_ + 
                             (1.0f - vocal_release_coeff_) * vocal_level;
        }
        
        // 楽器成分検出（主にSideチャンネル + Mid低域・高域）
        float inst_low = instrument_low_.process(mid);
        float inst_high = instrument_high_.process(mid);
        float inst_stereo = std::abs(side);
        float instrument_level = std::max({std::abs(inst_low), 
                                          std::abs(inst_high), 
                                          inst_stereo});
        
        if (instrument_level > instrument_envelope_) {
            instrument_envelope_ = inst_attack_coeff_ * instrument_envelope_ + 
                                  (1.0f - inst_attack_coeff_) * instrument_level;
        } else {
            instrument_envelope_ = inst_release_coeff_ * instrument_envelope_ + 
                                  (1.0f - inst_release_coeff_) * instrument_level;
        }
        
        // 動的分離処理
        return applyDynamicSeparation(mid, side);
    }
    
    std::pair<float, float> applyDynamicSeparation(float mid, float side) {
        // ボーカル優勢度の計算
        float vocal_dominance = vocal_envelope_ / (vocal_envelope_ + instrument_envelope_ + 1e-10f);
        float instrument_dominance = 1.0f - vocal_dominance;
        
        // 適応的処理
        float enhanced_mid = mid;
        float enhanced_side = side;
        
        // ボーカルが優勢な場合：Midを強化、Sideを抑制
        if (vocal_dominance > 0.6f) {
            float vocal_boost = 1.0f + vocal_enhance_ * vocal_dominance;
            enhanced_mid *= vocal_boost;
            enhanced_side *= (1.0f - vocal_enhance_ * 0.3f);
        }
        // 楽器が優勢な場合：Sideを強化、ステレオ幅拡張
        else if (instrument_dominance > 0.7f) {
            float stereo_boost = 1.0f + instrument_enhance_ * instrument_dominance;
            enhanced_side *= stereo_boost * stereo_width_;
        }
        
        return {enhanced_mid, enhanced_side};
    }
};

// 動的マルチバンドセパレーター
class DynamicMultibandSeparator {
public:
    struct Band {
        double freq_low, freq_high;
        bool vocal_dominant;  // この帯域はボーカル優勢か
        double separation_strength;
        SimpleBiquad lpf, hpf;
        float envelope = 0.0f;
        double attack_coeff = 0.0f;
        double release_coeff = 0.0f;
    };
    
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        enabled_ = params.value("enabled", true);
        
        // 帯域設定
        bands_ = {
            {60.0, 250.0, false, 0.1, {}, {}, 0.0f, 0.0, 0.0},      // 低域楽器
            {250.0, 800.0, false, 0.15, {}, {}, 0.0f, 0.0, 0.0},    // 中低域楽器
            {800.0, 2000.0, true, 0.3, {}, {}, 0.0f, 0.0, 0.0},     // ボーカル基音
            {2000.0, 5000.0, true, 0.4, {}, {}, 0.0f, 0.0, 0.0},    // ボーカル主要域
            {5000.0, 8000.0, false, 0.2, {}, {}, 0.0f, 0.0, 0.0},   // 楽器高域
            {8000.0, 20000.0, false, 0.1, {}, {}, 0.0f, 0.0, 0.0}   // 超高域
        };
        
        // 各帯域のフィルター初期化
        for (auto& band : bands_) {
            band.lpf.set_lpf(sr, band.freq_high, 0.707);
            band.hpf.set_hpf(sr, band.freq_low, 0.707);
            
            // エンベロープ特性
            double attack_ms = band.vocal_dominant ? 2.0 : 5.0;
            double release_ms = band.vocal_dominant ? 50.0 : 100.0;
            band.attack_coeff = std::exp(-1.0 / (attack_ms / 1000.0 * sr));
            band.release_coeff = std::exp(-1.0 / (release_ms / 1000.0 * sr));
        }
    }
    
    std::pair<float, float> process(float left, float right) {
        if (!enabled_) return {left, right};
        
        float enhanced_left = 0.0f;
        float enhanced_right = 0.0f;
        
        // 各帯域を処理
        for (auto& band : bands_) {
            auto band_signals = extractBand(left, right, band);
            auto processed_band = processBand(band_signals.first, band_signals.second, band);
            
            enhanced_left += processed_band.first;
            enhanced_right += processed_band.second;
        }
        
        return {enhanced_left, enhanced_right};
    }

private:
    double sample_rate_ = 44100.0;
    bool enabled_ = true;
    std::vector<Band> bands_;
    
    std::pair<float, float> extractBand(float left, float right, Band& band) {
        // 帯域抽出
        float band_left = band.hpf.process(left);
        band_left = band.lpf.process(band_left);
        
        float band_right = band.hpf.process(right);
        band_right = band.lpf.process(band_right);
        
        return {band_left, band_right};
    }
    
    std::pair<float, float> processBand(float left, float right, Band& band) {
        // M/S変換
        float mid = (left + right) * 0.5f;
        float side = (left - right) * 0.5f;
        
        // エンベロープ追跡
        float level = band.vocal_dominant ? std::abs(mid) : std::abs(side);
        
        if (level > band.envelope) {
            band.envelope = band.attack_coeff * band.envelope + 
                           (1.0f - band.attack_coeff) * level;
        } else {
            band.envelope = band.release_coeff * band.envelope + 
                           (1.0f - band.release_coeff) * level;
        }
        
        // 動的分離
        float separation_gain = 1.0f + band.separation_strength * band.envelope;
        
        if (band.vocal_dominant) {
            // ボーカル帯域：Midを強化
            mid *= separation_gain;
            side *= (2.0f - separation_gain) * 0.8f;  // Sideを少し抑制
        } else {
            // 楽器帯域：Sideを強化
            side *= separation_gain;
        }
        
        // L/R復元
        float processed_left = mid + side;
        float processed_right = mid - side;
        
        return {processed_left, processed_right};
    }
};

// スペクトラル・ボーカル分離（高度版）
class SpectralVocalSeparator {
public:
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        fft_size_ = params.value("fft_size", 2048);
        enabled_ = params.value("enabled", false);  // 計算負荷が高いためデフォルト無効
        
        vocal_threshold_ = params.value("vocal_threshold", 0.7);
        separation_strength_ = params.value("separation_strength", 0.3);
        
        // ウィンドウ関数
        window_.resize(fft_size_);
        for (size_t i = 0; i < fft_size_; ++i) {
            window_[i] = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (fft_size_ - 1)));
        }
        
        // バッファ初期化
        input_buffer_.resize(fft_size_, 0.0f);
        output_buffer_.resize(fft_size_, 0.0f);
        overlap_buffer_.resize(fft_size_, 0.0f);
    }
    
    std::pair<float, float> process(float left, float right) {
        if (!enabled_) return {left, right};
        
        // 簡易版実装（実際にはFFT処理が必要）
        // ここでは概念的な処理のみ示す
        
        // ボーカル検出（中央定位 + 特定周波数特性）
        float mid = (left + right) * 0.5f;
        float side = (left - right) * 0.5f;
        
        // 簡易ボーカル検出
        bool vocal_detected = detectVocalActivity(mid, side);
        
        if (vocal_detected) {
            // ボーカル強化モード
            mid *= (1.0f + separation_strength_);
            side *= (1.0f - separation_strength_ * 0.5f);
        } else {
            // 楽器強化モード
            side *= (1.0f + separation_strength_ * 0.7f);
        }
        
        return {mid + side, mid - side};
    }

private:
    double sample_rate_ = 44100.0;
    size_t fft_size_ = 2048;
    bool enabled_ = false;
    double vocal_threshold_ = 0.7;
    double separation_strength_ = 0.3;
    
    std::vector<float> window_;
    std::vector<float> input_buffer_, output_buffer_, overlap_buffer_;
    
    bool detectVocalActivity(float mid, float side) {
        // 簡易ボーカル検出アルゴリズム
        float center_energy = std::abs(mid);
        float stereo_energy = std::abs(side);
        
        // ボーカルは通常中央に定位
        float centrality = center_energy / (center_energy + stereo_energy + 1e-10f);
        
        return centrality > vocal_threshold_;
    }
};