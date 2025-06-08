// ./enhanced_separation_processor.h
#pragma once
#include "vocal_instrument_separator.h"

class EnhancedSeparationProcessor {
public:
    void setup(double sr, int channels, const json& params) {
        sample_rate_ = sr;
        channels_ = channels;
        
        if (channels != 2) {
            std::cout << "Warning: Separation requires stereo input. Proceeding without separation." << std::endl;
            separation_enabled_ = false;
            return;
        }
        
        separation_enabled_ = params.value("separation_enabled", true);
        
        if (separation_enabled_) {
            // 分離手法の選択
            std::string method = params.value("separation_method", "ms_adaptive");
            
            if (method == "ms_adaptive") {
                ms_separator_.setup(sr, params.value("ms_separator", json({})));
                use_ms_separator_ = true;
            } else if (method == "multiband") {
                multiband_separator_.setup(sr, params.value("multiband_separator", json({})));
                use_multiband_separator_ = true;
            } else if (method == "spectral") {
                spectral_separator_.setup(sr, params.value("spectral_separator", json({})));
                use_spectral_separator_ = true;
            } else if (method == "hybrid") {
                // 複数手法の組み合わせ
                ms_separator_.setup(sr, params.value("ms_separator", json({})));
                multiband_separator_.setup(sr, params.value("multiband_separator", json({})));
                use_ms_separator_ = true;
                use_multiband_separator_ = true;
            }
        }
        
        // 従来のエフェクト処理も継続
        setupTraditionalProcessors(sr, params);
    }
    
    void processBuffer(std::vector<float>& buffer) {
        if (channels_ != 2) {
            // モノラルまたは非対応チャンネル数
            processTraditionalEffects(buffer);
            return;
        }
        
        size_t frame_count = buffer.size() / 2;
        
        for (size_t i = 0; i < frame_count; ++i) {
            float left = buffer[i * 2];
            float right = buffer[i * 2 + 1];
            
            // 分離処理
            if (separation_enabled_) {
                auto separated = applySeparation(left, right);
                left = separated.first;
                right = separated.second;
            }
            
            // 従来のエフェクト処理
            auto processed = applyTraditionalEffects(left, right);
            buffer[i * 2] = processed.first;
            buffer[i * 2 + 1] = processed.second;
        }
    }

private:
    double sample_rate_ = 44100.0;
    int channels_ = 2;
    bool separation_enabled_ = true;
    
    // 分離プロセッサー
    MSVocalInstrumentSeparator ms_separator_;
    DynamicMultibandSeparator multiband_separator_;
    SpectralVocalSeparator spectral_separator_;
    
    bool use_ms_separator_ = false;
    bool use_multiband_separator_ = false;
    bool use_spectral_separator_ = false;
    
    // 従来のプロセッサー（モノラル処理用）
    SimpleExciterProcessor exciter_l_, exciter_r_;
    MultiParametricEQProcessor eq_processor_;
    SimpleLimiter final_limiter_;
    
    void setupTraditionalProcessors(double sr, const json& params) {
        exciter_l_.setup(sr, params.value("exciter", json({})));
        exciter_r_.setup(sr, params.value("exciter", json({})));
        eq_processor_.setup(sr, params.value("peq", json({})));
        final_limiter_.setup(sr);
    }
    
    std::pair<float, float> applySeparation(float left, float right) {
        std::pair<float, float> result = {left, right};
        
        // M/S分離適用
        if (use_ms_separator_) {
            result = ms_separator_.process(result.first, result.second);
        }
        
        // マルチバンド分離適用
        if (use_multiband_separator_) {
            result = multiband_separator_.process(result.first, result.second);
        }
        
        // スペクトラル分離適用
        if (use_spectral_separator_) {
            result = spectral_separator_.process(result.first, result.second);
        }
        
        return result;
    }
    
    std::pair<float, float> applyTraditionalEffects(float left, float right) {
        // チャンネルごとにエキサイター適用
        left = exciter_l_.process(left);
        right = exciter_r_.process(right);
        
        // ステレオ対応EQ処理
        std::vector<float> stereo_buffer = {left, right};
        eq_processor_.process_inplace(stereo_buffer, 2);
        
        return {stereo_buffer[0], stereo_buffer[1]};
    }
    
    void processTraditionalEffects(std::vector<float>& buffer) {
        // モノラル処理のフォールバック
        size_t frame_count = buffer.size() / channels_;
        
        for (size_t i = 0; i < frame_count; ++i) {
            for (int ch = 0; ch < channels_; ++ch) {
                float sample = buffer[i * channels_ + ch];
                if (ch == 0) {
                    sample = exciter_l_.process(sample);
                } else {
                    sample = exciter_r_.process(sample);
                }
                buffer[i * channels_ + ch] = sample;
            }
        }
        
        eq_processor_.process_inplace(buffer, channels_);
        final_limiter_.process_inplace(buffer, channels_);
    }
};

// 分離品質分析ツール
class SeparationAnalyzer {
public:
    struct SeparationMetrics {
        double vocal_clarity_score;
        double instrument_separation_score;
        double stereo_width_before;
        double stereo_width_after;
        double vocal_prominence_before;
        double vocal_prominence_after;
    };
    
    static SeparationMetrics analyze(const std::vector<float>& before_buffer,
                                    const std::vector<float>& after_buffer,
                                    int channels, double sample_rate) {
        SeparationMetrics metrics;
        
        if (channels != 2) {
            // ステレオでない場合は分析不可
            return metrics;
        }
        
        // 処理前後の分析
        auto before_analysis = analyzeBuffer(before_buffer, sample_rate);
        auto after_analysis = analyzeBuffer(after_buffer, sample_rate);
        
        // 分離効果の評価
        metrics.vocal_clarity_score = calculateVocalClarity(after_analysis);
        metrics.instrument_separation_score = calculateInstrumentSeparation(after_analysis);
        metrics.stereo_width_before = before_analysis.stereo_width;
        metrics.stereo_width_after = after_analysis.stereo_width;
        metrics.vocal_prominence_before = before_analysis.vocal_prominence;
        metrics.vocal_prominence_after = after_analysis.vocal_prominence;
        
        return metrics;
    }

private:
    struct BufferAnalysis {
        double stereo_width;
        double vocal_prominence;
        std::vector<double> frequency_content;
    };
    
    static BufferAnalysis analyzeBuffer(const std::vector<float>& buffer, double sample_rate) {
        BufferAnalysis analysis;
        
        size_t frame_count = buffer.size() / 2;
        double mid_energy = 0.0, side_energy = 0.0;
        double vocal_band_energy = 0.0, total_energy = 0.0;
        
        for (size_t i = 0; i < frame_count; ++i) {
            float left = buffer[i * 2];
            float right = buffer[i * 2 + 1];
            float mid = (left + right) * 0.5f;
            float side = (left - right) * 0.5f;
            
            mid_energy += mid * mid;
            side_energy += side * side;
            total_energy += left * left + right * right;
            
            // 簡易ボーカル帯域エネルギー（1-5kHz）
            // 実際にはバンドパスフィルターが必要
            vocal_band_energy += mid * mid; // 簡略化
        }
        
        analysis.stereo_width = side_energy / (mid_energy + side_energy + 1e-10);
        analysis.vocal_prominence = vocal_band_energy / (total_energy + 1e-10);
        
        return analysis;
    }
    
    static double calculateVocalClarity(const BufferAnalysis& analysis) {
        // ボーカルの明瞭度スコア計算
        return std::min(1.0, analysis.vocal_prominence * 2.0);
    }
    
    static double calculateInstrumentSeparation(const BufferAnalysis& analysis) {
        // 楽器分離スコア計算
        return std::min(1.0, analysis.stereo_width * 1.5);
    }
};