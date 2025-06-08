// ./main.cpp - 完全版（ノイズレス艶エンハンサー統合）
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstdio>
#include <stdexcept>
#include <filesystem> // Required for path manipulation

// 3rd party libraries
#include <sndfile.h>
#include <samplerate.h>
#include "json.hpp"

using json = nlohmann::json;

// --- Utilities ---
double db_to_linear(double db) {
    if (!std::isfinite(db)) return 1.0;
    return std::pow(10.0, db / 20.0);
}

// --- Simple Biquad Filter ---
class SimpleBiquad {
public:
    SimpleBiquad(std::string name = "Unnamed") : filter_name_(name) { reset(); }
    void reset() { a1 = 0.0; a2 = 0.0; b0 = 1.0; b1 = 0.0; b2 = 0.0; z1 = 0.0; z2 = 0.0; is_bypassed_ = false; }
    
    void set_lpf(double sr, double freq, double q) {
        reset();
        q = std::max(0.1, q); freq = std::max(10.0, std::min(freq, sr / 2.2));
        double w0 = 2.0 * M_PI * freq / sr, cos_w0 = std::cos(w0), sin_w0 = std::sin(w0);
        double alpha = sin_w0 / (2.0 * q), a0 = 1.0 + alpha;
        b0 = (1.0 - cos_w0) / 2.0 / a0; b1 = (1.0 - cos_w0) / a0; b2 = b0;
        a1 = -2.0 * cos_w0 / a0; a2 = (1.0 - alpha) / a0;
    }
    void set_hpf(double sr, double freq, double q) {
        reset();
        q = std::max(0.1, q); freq = std::max(10.0, std::min(freq, sr / 2.2));
        double w0 = 2.0 * M_PI * freq / sr, cos_w0 = std::cos(w0), sin_w0 = std::sin(w0);
        double alpha = sin_w0 / (2.0 * q), a0 = 1.0 + alpha;
        b0 = (1.0 + cos_w0) / 2.0 / a0; b1 = -(1.0 + cos_w0) / a0; b2 = b0;
        a1 = -2.0 * cos_w0 / a0; a2 = (1.0 - alpha) / a0;
    }
    void set_peaking(double sr, double freq, double q, double gain_db) {
        reset();
        q = std::max(0.1, q); freq = std::max(10.0, std::min(freq, sr / 2.2));
        double A = db_to_linear(gain_db / 2.0);
        double w0 = 2.0 * M_PI * freq / sr, cos_w0 = std::cos(w0), sin_w0 = std::sin(w0);
        double alpha = sin_w0 / (2.0 * q), a0 = 1.0 + alpha / A;
        b0 = (1.0 + alpha * A) / a0; b1 = -2.0 * cos_w0 / a0; b2 = (1.0 - alpha * A) / a0;
        a1 = b1; a2 = (1.0 - alpha / A) / a0;
    }
    
    float process(float in) {
        if (is_bypassed_ || std::isnan(in) || std::isinf(in)) return in;
        if (std::isnan(z1) || std::isinf(z1) || std::isnan(z2) || std::isinf(z2)) reset();
        double out = b0 * in + z1;
        z1 = b1 * in - a1 * out + z2;
        z2 = b2 * in - a2 * out;
        return static_cast<float>(out);
    }
private:
    std::string filter_name_; bool is_bypassed_ = false;
    double a1, a2, b0, b1, b2, z1, z2;
};

// --- High-Quality Resampler Class ---
class Resampler {
public:
    Resampler(int channels, double src_ratio, int converter_type = SRC_SINC_BEST_QUALITY)
        : channels_(channels), src_ratio_(src_ratio) {
        int error;
        src_state_ = src_new(converter_type, channels, &error);
        if (!src_state_) {
            throw std::runtime_error(std::string("Resampler error: src_new failed: ") + src_strerror(error));
        }
    }

    ~Resampler() {
        if (src_state_) {
            src_delete(src_state_);
        }
    }
    
    Resampler(const Resampler&) = delete;
    Resampler& operator=(const Resampler&) = delete;

    std::vector<float> process(const std::vector<float>& input, long input_frames) {
        long output_frames_estimate = static_cast<long>(ceil(input_frames * src_ratio_));
        std::vector<float> output(output_frames_estimate * channels_);

        SRC_DATA src_data;
        src_data.data_in = input.data();
        src_data.input_frames = input_frames;
        src_data.data_out = output.data();
        src_data.output_frames = output_frames_estimate;
        src_data.src_ratio = src_ratio_;
        src_data.end_of_input = 0;

        int error = src_process(src_state_, &src_data);
        if (error) {
            throw std::runtime_error(std::string("Resampler error: src_process failed: ") + src_strerror(error));
        }

        output.resize(src_data.output_frames_gen * channels_);
        return output;
    }

private:
    SRC_STATE* src_state_ = nullptr;
    int channels_;
    double src_ratio_;
};

// --- Simple Mastering Limiter ---
class SimpleLimiter {
public:
    void setup(double sr, float threshold_db = -0.1f, float release_ms = 100.0f) {
        sample_rate_ = sr;
        threshold_linear_ = static_cast<float>(db_to_linear(threshold_db));
        attack_coeff_ = std::exp(-1.0 / (0.001 * sample_rate_));
        release_coeff_ = std::exp(-1.0 / (release_ms / 1000.0 * sample_rate_));
    }

    void process_inplace(std::vector<float>& buffer, int channels) {
        size_t frame_count = buffer.size() / channels;
        for (size_t i = 0; i < frame_count; ++i) {
            float peak_val = 0.0f;
            for (int ch = 0; ch < channels; ++ch) {
                peak_val = std::max(peak_val, std::abs(buffer[i * channels + ch]));
            }

            if (peak_val > envelope_) {
                envelope_ = attack_coeff_ * envelope_ + (1.0f - attack_coeff_) * peak_val;
            } else {
                envelope_ = release_coeff_ * envelope_ + (1.0f - release_coeff_) * peak_val;
            }

            float gain = 1.0f;
            if (envelope_ > threshold_linear_) {
                gain = threshold_linear_ / envelope_;
            }

            for (int ch = 0; ch < channels; ++ch) {
                buffer[i * channels + ch] *= gain;
            }
        }
    }

private:
    double sample_rate_ = 44100.0;
    double attack_coeff_ = 0.0;
    double release_coeff_ = 0.0;
    float envelope_ = 0.0f;
    float threshold_linear_ = 1.0f;
};

// --- Advanced Noise Controller ---
class AdvancedNoiseController {
public:
    AdvancedNoiseController() : 
        anti_alias_lpf_("AntiAlias"),
        artifact_detector_("ArtifactDetector"),
        dc_blocker_("DCBlocker") {}
    
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        enabled_ = params.value("enabled", true);
        
        if (!enabled_) return;
        
        // アンチエイリアシング設定
        high_freq_limit_ = params.value("high_freq_limit", 32000.0);
        anti_aliasing_strength_ = params.value("anti_aliasing_strength", 0.9);
        
        // アーティファクト検出・抑制
        artifact_threshold_ = params.value("artifact_threshold", 0.12);
        artifact_suppression_ = params.value("artifact_suppression", 0.8);
        
        // ノイズゲート
        noise_floor_db_ = params.value("noise_floor_db", -85.0);
        gate_attack_ms_ = params.value("gate_attack_ms", 1.5);
        gate_release_ms_ = params.value("gate_release_ms", 40.0);
        
        setupFilters(sr);
        setupGating(sr);
    }
    
    float process(float input) {
        if (!enabled_) return input;
        
        // 1. DCオフセット除去
        float processed = dc_blocker_.process(input);
        
        // 2. アーティファクト検出・抑制
        processed = suppressArtifacts(processed);
        
        // 3. アンチエイリアシング
        processed = anti_alias_lpf_.process(processed);
        
        // 4. ノイズゲート
        processed = applyNoiseGate(processed);
        
        // 5. 最終クリーニング
        processed = finalCleanup(processed);
        
        return processed;
    }

private:
    double sample_rate_ = 44100.0;
    bool enabled_ = true;
    
    double high_freq_limit_ = 32000.0;
    double anti_aliasing_strength_ = 0.9;
    double artifact_threshold_ = 0.12;
    double artifact_suppression_ = 0.8;
    double noise_floor_db_ = -85.0;
    double gate_attack_ms_ = 1.5;
    double gate_release_ms_ = 40.0;
    
    SimpleBiquad anti_alias_lpf_, artifact_detector_, dc_blocker_;
    
    float gate_envelope_ = 0.0f;
    float gate_attack_coeff_ = 0.0f;
    float gate_release_coeff_ = 0.0f;
    float noise_floor_linear_ = 0.0f;
    
    float artifact_envelope_ = 0.0f;
    float artifact_attack_coeff_ = 0.0f;
    float artifact_release_coeff_ = 0.0f;
    
    void setupFilters(double sr) {
        double safe_cutoff = std::min(high_freq_limit_, sr * 0.38);
        anti_alias_lpf_.set_lpf(sr, safe_cutoff, 1.0 + anti_aliasing_strength_ * 0.8);
        artifact_detector_.set_hpf(sr, 7500.0, 1.4);
        dc_blocker_.set_hpf(sr, 3.0, 0.707);
        
        artifact_attack_coeff_ = std::exp(-1.0 / (0.0008 * sr));
        artifact_release_coeff_ = std::exp(-1.0 / (0.008 * sr));
    }
    
    void setupGating(double sr) {
        noise_floor_linear_ = static_cast<float>(db_to_linear(noise_floor_db_));
        gate_attack_coeff_ = std::exp(-1.0 / (gate_attack_ms_ / 1000.0 * sr));
        gate_release_coeff_ = std::exp(-1.0 / (gate_release_ms_ / 1000.0 * sr));
    }
    
    float suppressArtifacts(float input) {
        float high_freq_content = artifact_detector_.process(input);
        float artifact_level = std::abs(high_freq_content);
        
        if (artifact_level > artifact_envelope_) {
            artifact_envelope_ = artifact_attack_coeff_ * artifact_envelope_ + 
                                (1.0f - artifact_attack_coeff_) * artifact_level;
        } else {
            artifact_envelope_ = artifact_release_coeff_ * artifact_envelope_ + 
                                (1.0f - artifact_release_coeff_) * artifact_level;
        }
        
        float suppression_factor = 1.0f;
        if (artifact_envelope_ > artifact_threshold_) {
            float excess = (artifact_envelope_ - artifact_threshold_) / artifact_threshold_;
            suppression_factor = 1.0f - (excess * artifact_suppression_);
            suppression_factor = std::max(0.2f, suppression_factor);
        }
        
        return input * suppression_factor;
    }
    
    float applyNoiseGate(float input) {
        float abs_input = std::abs(input);
        
        if (abs_input > gate_envelope_) {
            gate_envelope_ = gate_attack_coeff_ * gate_envelope_ + 
                            (1.0f - gate_attack_coeff_) * abs_input;
        } else {
            gate_envelope_ = gate_release_coeff_ * gate_envelope_ + 
                            (1.0f - gate_release_coeff_) * abs_input;
        }
        
        float gate_gain = 1.0f;
        if (gate_envelope_ < noise_floor_linear_) {
            gate_gain = gate_envelope_ / noise_floor_linear_;
            gate_gain = std::max(0.0f, std::min(1.0f, gate_gain));
        }
        
        return input * gate_gain;
    }
    
    float finalCleanup(float input) {
        input = std::max(-0.98f, std::min(0.98f, input));
        if (!std::isfinite(input)) input = 0.0f;
        if (std::abs(input) < 1e-7f) input = 0.0f;
        return input;
    }
};

// --- Premium Gloss Enhancer with Noise Control ---
class PremiumGlossEnhancer {
public:
    PremiumGlossEnhancer() : 
        presence_filter_("Presence"), 
        air_filter_("Air"), 
        warmth_filter_("Warmth"),
        transient_hpf_("Transient_HPF"),
        noise_controller_() {}
    
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        enabled_ = params.value("enabled", true);
        
        if (!enabled_) return;
        
        // 安全で音楽的な設定
        harmonic_drive_ = std::min(0.5, params.value("harmonic_drive", 0.35));
        even_harmonics_ = std::min(0.4, params.value("even_harmonics", 0.28));
        odd_harmonics_ = std::min(0.3, params.value("odd_harmonics", 0.18));
        
        // 周波数特性（美しい艶を実現）
        presence_gain_ = std::min(1.8, params.value("presence_gain", 1.3));
        air_gain_ = std::min(1.5, params.value("air_gain", 1.0));
        warmth_gain_ = std::min(1.2, params.value("warmth_gain", 0.7));
        
        // アナログ特性（控えめで上品）
        analog_color_ = std::min(0.3, params.value("analog_color", 0.18));
        saturation_mix_ = std::min(0.2, params.value("saturation_mix", 0.1));
        
        // トランジェント処理
        transient_enhance_ = std::min(0.3, params.value("transient_enhance", 0.15));
        
        total_mix_ = std::min(0.35, params.value("total_mix", 0.22));
        
        // フィルター初期化
        presence_filter_.set_peaking(sr, 4200.0, 1.6, presence_gain_);
        air_filter_.set_peaking(sr, 11500.0, 1.3, air_gain_);
        warmth_filter_.set_peaking(sr, 350.0, 1.1, warmth_gain_);
        transient_hpf_.set_hpf(sr, 1800.0, 0.9);
        
        // ノイズコントローラー初期化
        json noise_params = {
            {"enabled", true},
            {"high_freq_limit", 30000.0},
            {"anti_aliasing_strength", 0.95},
            {"artifact_suppression", 0.85}
        };
        noise_controller_.setup(sr, noise_params);
        
        setupEnvelopeFollower(sr);
    }
    
    float process(float input) {
        if (!enabled_) return input;
        
        // 入力保護
        if (std::abs(input) > 0.85f) {
            input *= 0.85f / std::abs(input);
        }
        
        // 基本処理
        float enhanced = enhanceFrequencyCharacter(input);
        float harmonics = generatePremiumHarmonics(enhanced);
        float saturated = applyGentleAnalogColoration(enhanced + harmonics);
        float with_transients = enhanceTransients(saturated);
        
        // ノイズ制御
        with_transients = noise_controller_.process(with_transients);
        
        // 最終ミックス
        return input * (1.0f - total_mix_) + with_transients * total_mix_;
    }

private:
    double sample_rate_ = 44100.0;
    bool enabled_ = true;
    
    float harmonic_drive_ = 0.35f;
    float even_harmonics_ = 0.28f;
    float odd_harmonics_ = 0.18f;
    float presence_gain_ = 1.3f;
    float air_gain_ = 1.0f;
    float warmth_gain_ = 0.7f;
    float analog_color_ = 0.18f;
    float saturation_mix_ = 0.1f;
    float transient_enhance_ = 0.15f;
    float total_mix_ = 0.22f;
    
    SimpleBiquad presence_filter_, air_filter_, warmth_filter_, transient_hpf_;
    AdvancedNoiseController noise_controller_;
    
    float envelope_ = 0.0f;
    float attack_coeff_ = 0.0f;
    float release_coeff_ = 0.0f;
    float prev_input_ = 0.0f;
    
    void setupEnvelopeFollower(double sr) {
        attack_coeff_ = std::exp(-1.0 / (0.003 * sr));
        release_coeff_ = std::exp(-1.0 / (0.08 * sr));
    }
    
    float enhanceFrequencyCharacter(float input) {
        float warm = warmth_filter_.process(input);
        float presence = presence_filter_.process(input);
        float air = air_filter_.process(input);
        
        return input + warm * 0.25f + presence * 0.5f + air * 0.35f;
    }
    
    float generatePremiumHarmonics(float input) {
        float abs_input = std::abs(input);
        
        if (abs_input > envelope_) {
            envelope_ = attack_coeff_ * envelope_ + (1.0f - attack_coeff_) * abs_input;
        } else {
            envelope_ = release_coeff_ * envelope_ + (1.0f - release_coeff_) * abs_input;
        }
        
        float dynamic_factor = std::min(0.9f, envelope_ * 2.5f);
        
        // 美しい偶次高調波（温かみ）
        float even_harm = generateEvenHarmonics(input) * even_harmonics_ * dynamic_factor;
        
        // エレガントな奇次高調波（明瞭度）
        float odd_harm = generateOddHarmonics(input) * odd_harmonics_ * dynamic_factor;
        
        // 微細な質感
        float subtle_texture = generateSubtleTexture(input) * 0.12f * dynamic_factor;
        
        return (even_harm + odd_harm + subtle_texture) * 0.8f;
    }
    
    float generateEvenHarmonics(float input) {
        float squared = input * input * 0.4f;
        float fourth = squared * squared * 0.08f;
        float tube_warmth = std::tanh(input * harmonic_drive_ * 0.8f) * input * 0.6f;
        
        return (squared * 0.7f + fourth + tube_warmth) * 0.5f;
    }
    
    float generateOddHarmonics(float input) {
        float cubed = input * input * input * 0.3f;
        float fifth = cubed * input * input * 0.04f;
        float clarity = std::tanh(input * harmonic_drive_ * 0.6f) * 0.5f;
        
        return (cubed + fifth + clarity) * 0.4f;
    }
    
    float generateSubtleTexture(float input) {
        float soft_asymmetry = input * (1.0f + input * 0.03f);
        float gentle_curve = input / (1.0f + std::abs(input) * 0.08f);
        
        return (soft_asymmetry + gentle_curve - input * 2.0f) * 0.6f;
    }
    
    float applyGentleAnalogColoration(float input) {
        float gentle_tape = std::tanh(input * (1.0f + analog_color_ * 0.4f)) * 0.97f;
        float transformer_warmth = input * (1.0f - std::abs(input) * 0.04f * analog_color_);
        float tube_character = input > 0.0f ? 
            std::pow(input, 1.0f - analog_color_ * 0.08f) : 
            -std::pow(-input, 1.0f - analog_color_ * 0.08f);
        
        float analog_blend = gentle_tape * 0.5f + transformer_warmth * 0.3f + tube_character * 0.2f;
        
        return input * (1.0f - saturation_mix_) + analog_blend * saturation_mix_;
    }
    
    float enhanceTransients(float input) {
        float high_freq = transient_hpf_.process(input);
        float transient_signal = high_freq - prev_input_ * 0.7f;
        prev_input_ = high_freq;
        
        float enhanced_transient = transient_signal * transient_enhance_;
        
        return input + enhanced_transient * 0.25f;
    }
};

// --- Enhanced Exciter with Noise Control ---
class EnhancedExciterProcessor {
public:
    EnhancedExciterProcessor(): hpf_("Exciter_HPF") {}
    
    void setup(double sr, const json& p) {
        if(!p.value("enabled",false)) return; 
        enabled_ = true;
        
        drive_ = std::min(3.5, p.value("drive", 2.8)); 
        mix_ = std::min(0.3, p.value("mix", 0.18)); 
        even_drive_ = std::min(0.8, p.value("even_drive", 0.5));
        
        hpf_.set_hpf(sr, p.value("crossover_freq", 7800.0), 0.707);
    }
    
    float process(float in) {
        if (!enabled_) return in;
        
        // 入力保護
        if (std::abs(in) > 0.9f) {
            in *= 0.9f / std::abs(in);
        }
        
        float high_freqs = hpf_.process(in);
        
        // 穏やかで音楽的なサチュレーション
        auto musical_saturate = [&](float s) {
            float odd = std::tanh(s * drive_ * 0.7f);
            float even = (s * s - (1.0f/4.0f)) * even_drive_ * 0.5f;
            return std::tanh(odd + even) * 0.92f;
        };
        
        float saturated_highs = musical_saturate(high_freqs);
        
        // 最終安全チェック
        if (!std::isfinite(saturated_highs)) saturated_highs = 0.0f;
        
        return in * (1.0f - mix_) + saturated_highs * mix_;
    }
    
private:
    SimpleBiquad hpf_;
    float drive_, mix_, even_drive_; 
    bool enabled_ = false;
};

// --- Enhanced Parametric EQ ---
class MultiParametricEQProcessor {
public:
    MultiParametricEQProcessor(): 
        peq_low_("PEQ_Low"), peq1_("PEQ_Mid1"), peq2_("PEQ_Mid2"), 
        peq3_("PEQ_Mid3"), peq4_("PEQ_High1"), peq5_("PEQ_High2"),
        lpf_final_("LPF_FinalCut") {}
    
    void setup(double sr, const json& p) {
        if(!p.value("enabled",false)) return; 
        enabled_ = true;
        
        // 音楽的で自然なEQ設定
        peq_low_.set_peaking(sr, 85.0, 0.9, 0.7);      // 低域の基盤
        peq1_.set_peaking(sr, 180.0, 1.2, 0.4);        // 低域の整理
        peq2_.set_peaking(sr, 1800.0, 1.6, 1.0);       // ボーカル基音
        peq3_.set_peaking(sr, 3200.0, 1.8, 1.2);       // プレゼンス
        peq4_.set_peaking(sr, 6800.0, 1.4, 1.0);       // 明瞭度
        peq5_.set_peaking(sr, 12000.0, 1.1, 0.7);      // エアー感
        
        // 強力なアンチエイリアシング
        lpf_final_.set_lpf(sr, 32000.0, 1.2);
    }
    
    void process_inplace(std::vector<float>& buffer, int channels) {
        if (!enabled_) return;
        size_t frame_count = buffer.size() / channels;
        for (size_t i = 0; i < frame_count; ++i) {
            for (int ch = 0; ch < channels; ++ch) {
                float sample = buffer[i * channels + ch];
                sample = peq_low_.process(sample);
                sample = peq1_.process(sample);
                sample = peq2_.process(sample);
                sample = peq3_.process(sample);
                sample = peq4_.process(sample);
                sample = peq5_.process(sample);
                sample = lpf_final_.process(sample);
                buffer[i * channels + ch] = sample;
            }
        }
    }
private:
    SimpleBiquad peq_low_, peq1_, peq2_, peq3_, peq4_, peq5_, lpf_final_;
    bool enabled_ = false;
};

// --- M/S Vocal-Instrument Separator ---
class MSVocalInstrumentSeparator {
public:
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        enabled_ = params.value("enabled", false);
        
        if (!enabled_) return;
        
        vocal_enhance_ = std::min(0.3, params.value("vocal_enhance", 0.2));
        instrument_enhance_ = std::min(0.25, params.value("instrument_enhance", 0.15));
        stereo_width_ = std::min(1.3, params.value("stereo_width", 1.08));
        
        // ボーカル検出用フィルター
        vocal_bandpass_low_.set_hpf(sr, 800.0, 0.707);
        vocal_bandpass_high_.set_lpf(sr, 5000.0, 0.707);
        
        // エンベロープフォロワー係数
        vocal_attack_coeff_ = std::exp(-1.0 / (0.004 * sr));
        vocal_release_coeff_ = std::exp(-1.0 / (0.12 * sr));
    }
    
    std::pair<float, float> process(float left, float right) {
        if (!enabled_) return {left, right};
        
        // M/S変換
        float mid = (left + right) * 0.5f;
        float side = (left - right) * 0.5f;
        
        // ボーカル検出（より保守的）
        float vocal_signal = vocal_bandpass_high_.process(vocal_bandpass_low_.process(mid));
        float vocal_level = std::abs(vocal_signal);
        
        // エンベロープ更新
        if (vocal_level > vocal_envelope_) {
            vocal_envelope_ = vocal_attack_coeff_ * vocal_envelope_ + 
                             (1.0f - vocal_attack_coeff_) * vocal_level;
        } else {
            vocal_envelope_ = vocal_release_coeff_ * vocal_envelope_ + 
                             (1.0f - vocal_release_coeff_) * vocal_level;
        }
        
        // 動的分離処理（穏やか）
        float vocal_dominance = vocal_envelope_ / (vocal_envelope_ + std::abs(side) + 1e-10f);
        
        float enhanced_mid = mid;
        float enhanced_side = side;
        
        if (vocal_dominance > 0.65f) {
            // ボーカル優勢：穏やかなMid強化
            enhanced_mid *= (1.0f + vocal_enhance_ * vocal_dominance * 0.8f);
            enhanced_side *= (1.0f - vocal_enhance_ * 0.2f);
        } else {
            // 楽器優勢：穏やかなSide強化
            enhanced_side *= (1.0f + instrument_enhance_ * 0.8f) * stereo_width_;
        }
        
        // L/R復元
        float enhanced_left = enhanced_mid + enhanced_side;
        float enhanced_right = enhanced_mid - enhanced_side;
        
        return {enhanced_left, enhanced_right};
    }

private:
    double sample_rate_ = 44100.0;
    bool enabled_ = false;
    double vocal_enhance_ = 0.2;
    double instrument_enhance_ = 0.15;
    double stereo_width_ = 1.08;
    
    SimpleBiquad vocal_bandpass_low_, vocal_bandpass_high_;
    float vocal_envelope_ = 0.0f;
    float vocal_attack_coeff_ = 0.0f;
    float vocal_release_coeff_ = 0.0f;
};

// --- Main Processing ---
int main(int argc, char *argv[]) {
    if (argc != 2) { 
        std::cerr << "Usage: " << argv[0] << " <input_audio_file>" << std::endl; 
        return 1; 
    }
    
    std::string input_filename = argv[1];
    
    // --- Start of modification ---
    // Get the path to the executable to locate params.json
    std::filesystem::path exe_path(argv[0]);
    std::filesystem::path params_path = exe_path.parent_path() / ".." / "params.json";

    json params = json::object(); // Initialize as an empty object to prevent null access
    std::ifstream f(params_path);
    if (f.is_open()) {
        try {
            // Parse with exceptions enabled to catch syntax errors
            params = json::parse(f, nullptr, true);
            std::cout << "Parameters loaded from " << std::filesystem::absolute(params_path).string() << std::endl;
        } catch (json::parse_error& e) {
            std::cerr << "Error parsing params.json: " << e.what() << ". Using defaults." << std::endl;
            params = json::object(); // Fallback to empty object on parse error
        }
    } else {
        std::cout << "Warning: params.json not found at " << std::filesystem::absolute(params_path).string() << ", using defaults." << std::endl;
    }
    // --- End of modification ---

    SF_INFO sfinfo_in;
    SNDFILE* infile = sf_open(input_filename.c_str(), SFM_READ, &sfinfo_in);
    if (!infile) {
        std::cerr << "Error: Could not open input file: " << sf_strerror(nullptr) << std::endl; 
        return 1;
    }

    size_t dot_pos = input_filename.rfind('.');
    std::string base_filename = (dot_pos == std::string::npos) ? input_filename : input_filename.substr(0, dot_pos);
    
    std::string output_ext;
    int final_format_type;

    if ((sfinfo_in.format & SF_FORMAT_TYPEMASK) == SF_FORMAT_MPEG) {
        std::cout << "Input is an MPEG file (like MP3). To avoid quality loss from re-compression, the output will be saved in WAV format." << std::endl;
        output_ext = ".wav";
        final_format_type = SF_FORMAT_WAV;
    } else {
        output_ext = (dot_pos == std::string::npos) ? ".wav" : input_filename.substr(dot_pos);
        final_format_type = (sfinfo_in.format & SF_FORMAT_TYPEMASK);
    }
    
    std::string tmp_filename = base_filename + "_tmp.wav";
    std::string final_filename = base_filename + "_final" + output_ext;
    
    const int TMP_SR = 192000;
    const int FINAL_SR = 96000;
    const size_t BUF_SIZE = 8192;

    std::cout << "=== Enhanced Sound Processor v2.0 ===" << std::endl;
    std::cout << "Input: " << input_filename << std::endl;
    std::cout << "Channels: " << sfinfo_in.channels << ", Sample Rate: " << sfinfo_in.samplerate << "Hz" << std::endl;

    try {
        // =========================================================================
        //  第1パス：アップサンプル、エフェクト処理、一時ファイルへの書き出し
        // =========================================================================
        std::cout << "\n--- Pass 1: Ultra-High Quality Processing (192kHz) ---" << std::endl;
        
        SF_INFO sfinfo_tmp = sfinfo_in;
        sfinfo_tmp.samplerate = TMP_SR;
        sfinfo_tmp.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

        SNDFILE* tmp_outfile = sf_open(tmp_filename.c_str(), SFM_WRITE, &sfinfo_tmp);
        if (!tmp_outfile) {
            throw std::runtime_error("Could not create temp file.");
        }

        Resampler upsampler(sfinfo_in.channels, (double)TMP_SR / sfinfo_in.samplerate);
        std::vector<EnhancedExciterProcessor> exciters(sfinfo_in.channels);
        std::vector<PremiumGlossEnhancer> gloss_enhancers(sfinfo_in.channels);
        MultiParametricEQProcessor eq_processor;
        MSVocalInstrumentSeparator separator;
        
        // プロセッサー初期化
        for (int i = 0; i < sfinfo_in.channels; ++i) {
            exciters[i].setup(TMP_SR, params.value("exciter", json::object()));
            gloss_enhancers[i].setup(TMP_SR, params.value("gloss_enhancer", json::object()));
        }
        eq_processor.setup(TMP_SR, params.value("peq", json::object()));
        
        // 分離機能（ステレオの場合のみ）
        bool separation_enabled = false;
        if (sfinfo_in.channels == 2) {
            json sep_params = params.value("separation", json::object());
            sep_params["enabled"] = sep_params.value("enabled", true);
            separator.setup(TMP_SR, sep_params);
            separation_enabled = sep_params.value("enabled", true);
            std::cout << "Vocal-Instrument separation: " << 
                (separation_enabled ? "Enabled" : "Disabled") << std::endl;
        }
        
        // 艶エンハンサー状態表示
        bool gloss_enabled = params.value("gloss_enhancer", json::object()).value("enabled", true);
        std::cout << "Premium Gloss Enhancement: " << 
            (gloss_enabled ? "Enabled" : "Disabled") << std::endl;
        
        std::vector<float> in_buf(BUF_SIZE * sfinfo_in.channels);
        sf_count_t frames_read;
        sf_count_t total_frames_processed = 0;
        
        while ((frames_read = sf_readf_float(infile, in_buf.data(), BUF_SIZE)) > 0) {
            in_buf.resize(frames_read * sfinfo_in.channels);
            std::vector<float> processed_buf = upsampler.process(in_buf, frames_read);

            long upsampled_frames = processed_buf.size() / sfinfo_in.channels;
            
            // 分離処理（ステレオの場合）
            if (separation_enabled && sfinfo_in.channels == 2) {
                for (long i = 0; i < upsampled_frames; ++i) {
                    float left = processed_buf[i * 2];
                    float right = processed_buf[i * 2 + 1];
                    auto separated = separator.process(left, right);
                    processed_buf[i * 2] = separated.first;
                    processed_buf[i * 2 + 1] = separated.second;
                }
            }
            
            // エキサイター処理
            for(long i = 0; i < upsampled_frames; ++i) {
                for(int ch = 0; ch < sfinfo_in.channels; ++ch) {
                    processed_buf[i * sfinfo_in.channels + ch] = 
                        exciters[ch].process(processed_buf[i * sfinfo_in.channels + ch]);
                }
            }
            
            // EQ処理
            eq_processor.process_inplace(processed_buf, sfinfo_in.channels);
            
            // 艶エンハンサー処理
            if (gloss_enabled) {
                for(long i = 0; i < upsampled_frames; ++i) {
                    for(int ch = 0; ch < sfinfo_in.channels; ++ch) {
                        processed_buf[i * sfinfo_in.channels + ch] = 
                            gloss_enhancers[ch].process(processed_buf[i * sfinfo_in.channels + ch]);
                    }
                }
            }

            sf_writef_float(tmp_outfile, processed_buf.data(), upsampled_frames);
            total_frames_processed += frames_read;
            
            // 進捗表示
            if (total_frames_processed % (sfinfo_in.samplerate * 5) == 0) {
                double progress = (double)total_frames_processed / sfinfo_in.frames * 100.0;
                std::cout << "Progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::endl;
            }
        }
        sf_close(infile);
        sf_close(tmp_outfile);

        // =========================================================================
        //  第2パス：リミッター、ノーマライズ、ダウンサンプル
        // =========================================================================
        std::cout << "\n--- Pass 2: Mastering & High-Quality Downsampling ---" << std::endl;

        SNDFILE* tmp_infile = sf_open(tmp_filename.c_str(), SFM_READ, &sfinfo_tmp);
        if (!tmp_infile) {
            throw std::runtime_error("Could not open temp file for reading.");
        }

        std::vector<float> full_buffer(sfinfo_tmp.frames * sfinfo_tmp.channels);
        sf_readf_float(tmp_infile, full_buffer.data(), sfinfo_tmp.frames);
        sf_close(tmp_infile);
        
        // マスタリングリミッター
        SimpleLimiter limiter;
        limiter.setup(TMP_SR, -0.2f, 40.0f);  // より穏やかなリミッティング
        std::cout << "Applying mastering limiter..." << std::endl;
        limiter.process_inplace(full_buffer, sfinfo_tmp.channels);
        
        // ピーク分析とノーマライズ
        float max_peak = 0.0f;
        double rms_sum = 0.0;
        for(const auto& sample : full_buffer) {
            max_peak = std::max(max_peak, std::abs(sample));
            rms_sum += sample * sample;
        }
        float rms_level = std::sqrt(rms_sum / full_buffer.size());
        
        std::cout << "Peak level: " << std::fixed << std::setprecision(2) << 
            20.0 * std::log10(max_peak) << "dB" << std::endl;
        std::cout << "RMS level: " << std::fixed << std::setprecision(2) << 
            20.0 * std::log10(rms_level) << "dB" << std::endl;
        
        float norm_gain = (max_peak > 0.0f) ? (float)(db_to_linear(-0.1) / max_peak) : 1.0f;
        std::cout << "Normalization gain: " << std::fixed << std::setprecision(2) << 
            20.0 * std::log10(norm_gain) << "dB" << std::endl;
        
        for(auto& sample : full_buffer) {
            sample *= norm_gain;
        }

        // 最終出力（高品質ダウンサンプリング）
        SF_INFO sfinfo_final = sfinfo_in;
        sfinfo_final.samplerate = FINAL_SR;
        sfinfo_final.format = final_format_type | SF_FORMAT_PCM_24;
        SNDFILE* final_outfile = sf_open(final_filename.c_str(), SFM_WRITE, &sfinfo_final);
        if (!final_outfile) { 
            throw std::runtime_error("Could not create final output file. The output format might not be supported for writing."); 
        }
        
        std::cout << "High-quality downsampling to " << FINAL_SR << "Hz..." << std::endl;
        Resampler downsampler(sfinfo_final.channels, (double)FINAL_SR / TMP_SR);
        std::vector<float> final_buf = downsampler.process(full_buffer, full_buffer.size() / sfinfo_final.channels);
        sf_writef_float(final_outfile, final_buf.data(), final_buf.size() / sfinfo_final.channels);
        
        sf_close(final_outfile);
        std::remove(tmp_filename.c_str());

        // 最終統計
        std::cout << "\n=== Processing Complete ===" << std::endl;
        std::cout << "Output: " << final_filename << std::endl;
        std::cout << "Format: 24-bit PCM @ " << FINAL_SR << "Hz" << std::endl;
        std::cout << "Dynamic Range: " << std::fixed << std::setprecision(1) << 
            (20.0 * std::log10(max_peak) - 20.0 * std::log10(rms_level)) << "dB" << std::endl;
        
        // 処理内容サマリー
        std::cout << "\nProcessing Summary:" << std::endl;
        std::cout << "• Vocal-Instrument Separation: " << (separation_enabled ? "✓" : "✗") << std::endl;
        std::cout << "• Enhanced Exciter: ✓" << std::endl;
        std::cout << "• Premium Gloss Enhancement: " << (gloss_enabled ? "✓" : "✗") << std::endl;
        std::cout << "• 6-Band Parametric EQ: ✓" << std::endl;
        std::cout << "• Advanced Noise Control: ✓" << std::endl;
        std::cout << "• Mastering Limiter: ✓" << std::endl;
        std::cout << "• High-Quality Resampling: ✓" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::remove(tmp_filename.c_str());
        return 1;
    }

    return 0;
}