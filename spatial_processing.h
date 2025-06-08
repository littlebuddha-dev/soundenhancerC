// ./spatial_processing.h
#pragma once
#include <vector>
#include <cmath>
#include <array>

// ステレオ幅調整とイメージング
class StereoEnhancer {
public:
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        width_ = params.value("width", 1.2); // 1.0 = normal, >1.0 = wider
        bass_mono_freq_ = params.value("bass_mono_freq", 120.0);
        enabled_ = params.value("enabled", true);
        
        // 低域モノ化フィルター
        bass_lpf_l_.set_lpf(sr, bass_mono_freq_, 0.707);
        bass_lpf_r_.set_lpf(sr, bass_mono_freq_, 0.707);
        bass_hpf_l_.set_hpf(sr, bass_mono_freq_, 0.707);
        bass_hpf_r_.set_hpf(sr, bass_mono_freq_, 0.707);
    }
    
    std::pair<float, float> process(float left, float right) {
        if (!enabled_) return {left, right};
        
        // M/S変換
        float mid = (left + right) * 0.5f;
        float side = (left - right) * 0.5f;
        
        // 低域処理（モノ化）
        float bass_l = bass_lpf_l_.process(left);
        float bass_r = bass_lpf_r_.process(right);
        float bass_mono = (bass_l + bass_r) * 0.5f;
        
        float high_l = bass_hpf_l_.process(left);
        float high_r = bass_hpf_r_.process(right);
        
        // ステレオ幅調整（高域のみ）
        float high_mid = (high_l + high_r) * 0.5f;
        float high_side = (high_l - high_r) * 0.5f * width_;
        
        // 再構成
        float processed_l = bass_mono + high_mid + high_side;
        float processed_r = bass_mono + high_mid - high_side;
        
        return {processed_l, processed_r};
    }

private:
    double sample_rate_ = 44100.0;
    double width_ = 1.2;
    double bass_mono_freq_ = 120.0;
    bool enabled_ = true;
    SimpleBiquad bass_lpf_l_, bass_lpf_r_, bass_hpf_l_, bass_hpf_r_;
};

// リバーブ（簡易版）
class SimpleReverb {
public:
    void setup(double sr, const json& params) {
        sample_rate_ = sr;
        room_size_ = params.value("room_size", 0.3);
        damping_ = params.value("damping", 0.5);
        wet_level_ = params.value("wet_level", 0.15);
        enabled_ = params.value("enabled", false);
        
        // 遅延線初期化
        setupDelayLines();
        
        // フィードバックフィルター
        for (auto& filter : damping_filters_) {
            filter.set_lpf(sr, 8000.0 * (1.0 - damping_), 0.707);
        }
    }
    
    std::pair<float, float> process(float left, float right) {
        if (!enabled_) return {left, right};
        
        float input = (left + right) * 0.5f;
        
        // オールパスフィルター段
        float allpass_out = input;
        for (size_t i = 0; i < allpass_delays_.size(); ++i) {
            allpass_out = processAllpass(allpass_out, i);
        }
        
        // 並列コムフィルター
        float reverb_l = 0.0f, reverb_r = 0.0f;
        for (size_t i = 0; i < comb_delays_l_.size(); ++i) {
            reverb_l += processComb(allpass_out, i, true);
            reverb_r += processComb(allpass_out, i, false);
        }
        
        // ステレオ装飾
        reverb_l = stereo_decorr_l_.process(reverb_l);
        reverb_r = stereo_decorr_r_.process(reverb_r);
        
        // ミックス
        float out_l = left * (1.0f - wet_level_) + reverb_l * wet_level_;
        float out_r = right * (1.0f - wet_level_) + reverb_r * wet_level_;
        
        return {out_l, out_r};
    }

private:
    double sample_rate_ = 44100.0;
    double room_size_ = 0.3;
    double damping_ = 0.5;
    double wet_level_ = 0.15;
    bool enabled_ = false;
    
    // 遅延線
    std::vector<std::vector<float>> allpass_delays_;
    std::vector<std::vector<float>> comb_delays_l_, comb_delays_r_;
    std::vector<size_t> allpass_indices_, comb_indices_l_, comb_indices_r_;
    std::vector<SimpleBiquad> damping_filters_;
    SimpleBiquad stereo_decorr_l_, stereo_decorr_r_;
    
    void setupDelayLines() {
        // Freeverb風のパラメータ
        std::vector<int> allpass_lengths = {347, 113, 37, 59};
        std::vector<int> comb_lengths_l = {1557, 1617, 1491, 1422, 1277, 1356, 1188, 1116};
        std::vector<int> comb_lengths_r = {1617, 1557, 1422, 1491, 1356, 1277, 1116, 1188};
        
        // サンプルレートに応じてスケール
        double scale = sample_rate_ / 44100.0;
        
        allpass_delays_.resize(allpass_lengths.size());
        allpass_indices_.resize(allpass_lengths.size(), 0);
        for (size_t i = 0; i < allpass_lengths.size(); ++i) {
            int len = static_cast<int>(allpass_lengths[i] * scale);
            allpass_delays_[i].resize(len, 0.0f);
        }
        
        comb_delays_l_.resize(comb_lengths_l.size());
        comb_delays_r_.resize(comb_lengths_r.size());
        comb_indices_l_.resize(comb_lengths_l.size(), 0);
        comb_indices_r_.resize(comb_lengths_r.size(), 0);
        
        for (size_t i = 0; i < comb_lengths_l.size(); ++i) {
            int len_l = static_cast<int>(comb_lengths_l[i] * scale);
            int len_r = static_cast<int>(comb_lengths_r[i] * scale);
            comb_delays_l_[i].resize(len_l, 0.0f);
            comb_delays_r_[i].resize(len_r, 0.0f);
        }
        
        damping_filters_.resize(comb_lengths_l.size());
        
        // ステレオ装飾フィルター
        stereo_decorr_l_.set_peaking(sample_rate_, 3000.0, 0.8, 1.5);
        stereo_decorr_r_.set_peaking(sample_rate_, 4000.0, 0.9, -1.2);
    }
    
    float processAllpass(float input, size_t index) {
        size_t& pos = allpass_indices_[index];
        std::vector<float>& delay = allpass_delays_[index];
        
        float delay_out = delay[pos];
        float feedback = -input * 0.5f + delay_out;
        delay[pos] = input + feedback * 0.5f;
        
        pos = (pos + 1) % delay.size();
        return feedback;
    }
    
    float processComb(float input, size_t index, bool is_left) {
        auto& delays = is_left ? comb_delays_l_ : comb_delays_r_;
        auto& indices = is_left ? comb_indices_l_ : comb_indices_r_;
        
        size_t& pos = indices[index];
        std::vector<float>& delay = delays[index];
        
        float delay_out = delay[pos];
        float feedback = damping_filters_[index].process(delay_out * room_size_);
        delay[pos] = input + feedback;
        
        pos = (pos + 1) % delay.size();
        return delay_out;
    }
};