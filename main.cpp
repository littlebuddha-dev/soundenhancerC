// ./enhanced_main.cpp
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
<cmath>
#include <algorithm>
#include <iomanip>
#include <cstdio>
#include <stdexcept>
#include <memory>

// 3rd party libraries
#include <sndfile.h>
#include <samplerate.h>
#include "json.hpp"

// 拡張プロセッサーをインクルード
#include "advanced_dynamics.h"
#include "spatial_processing.h"
#include "advanced_eq_harmonics.h"

using json = nlohmann::json;

// --- 改良版オーディオプロセッサー ---
class EnhancedAudioProcessor {
public:
    void setup(double sr, int channels, const json& params) {
        sample_rate_ = sr;
        channels_ = channels;
        
        // 各プロセッサーの初期化
        spectral_gate_.setup(sr, params.value("spectral_gate", json({})));
        
        // マルチバンドコンプレッサーのセットアップ
        setupMultibandCompressor(params.value("multiband_comp", json({})));
        
        analog_saturation_.setup(sr, params.value("analog_saturation", json({})));
        harmonic_enhancer_.setup(sr, params.value("harmonic_enhancer", json({})));
        
        // ステレオプロセッサー（ステレオの場合のみ）
        if (channels == 2) {
            stereo_enhancer_.setup(sr, params.value("stereo_enhancer", json({})));
            reverb_.setup(sr, params.value("reverb", json({})));
        }
        
        // 最終段リミッター
        final_limiter_.setup(sr, 
            params.value("final_limiter", json({})).value("threshold_db", -0.3),
            params.value("final_limiter", json({})).value("release_ms", 50.0));
    }
    
    void processBuffer(std::vector<float>& buffer) {
        size_t frame_count = buffer.size() / channels_;
        
        for (size_t i = 0; i < frame_count; ++i) {
            if (channels_ == 1) {
                // モノラル処理
                float sample = buffer[i];
                sample = processMono(sample);
                buffer[i] = sample;
            } else if (channels_ == 2) {
                // ステレオ処理
                float left = buffer[i * 2];
                float right = buffer[i * 2 + 1];
                auto processed = processStereo(left, right);
                buffer[i * 2] = processed.first;
                buffer[i * 2 + 1] = processed.second;
            }
        }
        
        // 最終段処理
        final_limiter_.process_inplace(buffer, channels_);
    }

private:
    double sample_rate_ = 44100.0;
    int channels_ = 2;
    
    // プロセッサー群
    SpectralGate spectral_gate_;
    MultibandCompressor multiband_comp_;
    AnalogSaturation analog_saturation_;
    HarmonicEnhancer harmonic_enhancer_;
    StereoEnhancer stereo_enhancer_;
    SimpleReverb reverb_;
    SimpleLimiter final_limiter_;
    
    void setupMultibandCompressor(const json& params) {
        if (!params.value("enabled", false)) return;
        
        // 4バンド構成のデフォルト設定
        std::vector<MultibandCompressor::Band> bands = {
            // 低域 (20Hz - 250Hz)
            {20.0, 250.0, -12.0, 3.0, 10.0, 100.0, 2.0, true, 0.0, 0.0, 0.0, {}, {}, {}},
            // 低中域 (250Hz - 2kHz) 
            {250.0, 2000.0, -8.0, 2.5, 5.0, 80.0, 1.0, true, 0.0, 0.0, 0.0, {}, {}, {}},
            // 高中域 (2kHz - 8kHz)
            {2000.0, 8000.0, -6.0, 2.0, 3.0, 60.0, 1.5, true, 0.0, 0.0, 0.0, {}, {}, {}},
            // 高域 (8kHz - 20kHz)
            {8000.0, 20000.0, -4.0, 1.8, 2.0, 40.0, 2.0, true, 0.0, 0.0, 0.0, {}, {}, {}}
        };
        
        // パラメータから設定を上書き
        auto bands_json = params.value("bands", json::array());
        for (size_t i = 0; i < std::min(bands.size(), bands_json.size()); ++i) {
            auto& band = bands[i];
            auto& band_json = bands_json[i];
            
            band.threshold_db = band_json.value("threshold_db", band.threshold_db);
            band.ratio = band_json.value("ratio", band.ratio);
            band.attack_ms = band_json.value("attack_ms", band.attack_ms);
            band.release_ms = band_json.value("release_ms", band.release_ms);
            band.makeup_gain_db = band_json.value("makeup_gain_db", band.makeup_gain_db);
            band.enabled = band_json.value("enabled", band.enabled);
        }
        
        multiband_comp_.setup(sample_rate_, bands);
    }
    
    float processMono(float input) {
        float processed = input;
        
        // 1. スペクトラルゲート（ノイズ除去）
        processed = spectral_gate_.process(processed);
        
        // 2. マルチバンドコンプレッション
        processed = multiband_comp_.process(processed);
        
        // 3. アナログサチュレーション
        processed = analog_saturation_.process(processed);
        
        // 4. ハーモニックエンハンサー
        processed = harmonic_enhancer_.process(processed);
        
        return processed;
    }
    
    std::pair<float, float> processStereo(float left, float right) {
        // モノラル処理を各チャンネルに適用
        left = processMono(left);
        right = processMono(right);
        
        // ステレオ処理
        auto stereo_enhanced = stereo_enhancer_.process(left, right);
        auto reverb_processed = reverb_.process(stereo_enhanced.first, stereo_enhanced.second);
        
        return reverb_processed;
    }
};

// --- 音質分析ツール ---
class AudioAnalyzer {
public:
    struct AnalysisResult {
        double peak_db;
        double rms_db;
        double dynamic_range;
        double stereo_width;  // ステレオの場合のみ
        std::vector<double> frequency_response; // FFT結果
    };
    
    static AnalysisResult analyze(const std::vector<float>& buffer, int channels, double sample_rate) {
        AnalysisResult result;
        
        // ピーク値とRMS計算
        float peak = 0.0f, rms_sum = 0.0f;
        for (float sample : buffer) {
            peak = std::max(peak, std::abs(sample));
            rms_sum += sample * sample;
        }
        
        result.peak_db = 20.0 * std::log10(std::max(peak, 1e-10f));
        result.rms_db = 10.0 * std::log10(rms_sum / buffer.size());
        result.dynamic_range = result.peak_db - result.rms_db;
        
        // ステレオ幅計算（ステレオの場合）
        if (channels == 2) {
            double correlation = 0.0;
            size_t frame_count = buffer.size() / 2;
            
            for (size_t i = 0; i < frame_count; ++i) {
                correlation += buffer[i * 2] * buffer[i * 2 + 1];
            }
            correlation /= frame_count;
            result.stereo_width = 1.0 - correlation; // 簡易的な幅計算
        }
        
        return result;
    }
};

// --- 改良版メイン関数 ---
int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_audio_file> [preset_name]" << std::endl;
        std::cerr << "Available presets: master, vocal, instrumental, podcast" << std::endl;
        return 1;
    }
    
    std::string input_filename = argv[1];
    std::string preset_name = (argc > 2) ? argv[2] : "master";
    std::string base_filename = input_filename.substr(0, input_filename.rfind('.'));
    std::string output_filename = base_filename + "_enhanced.wav";
    
    // プリセット読み込み
    json params = loadPreset(preset_name);
    
    // 入力ファイルを開く
    SF_INFO sfinfo_in;
    SNDFILE* infile = sf_open(input_filename.c_str(), SFM_READ, &sfinfo_in);
    if (!infile) {
        std::cerr << "Error: Could not open input file." << std::endl;
        return 1;
    }
    
    std::cout << "Input: " << input_filename << std::endl;
    std::cout << "Channels: " << sfinfo_in.channels << ", Sample Rate: " << sfinfo_in.samplerate << "Hz" << std::endl;
    std::cout << "Preset: " << preset_name << std::endl;
    
    const int PROCESSING_SR = 192000;  // 高品質処理用サンプルレート
    const int OUTPUT_SR = 96000;       // 出力サンプルレート
    const size_t BUFFER_SIZE = 16384;  // 大きなバッファでより効率的に
    
    try {
        // アップサンプラー
        Resampler upsampler(sfinfo_in.channels, (double)PROCESSING_SR / sfinfo_in.samplerate);
        
        // 処理エンジン初期化
        EnhancedAudioProcessor processor;
        processor.setup(PROCESSING_SR, sfinfo_in.channels, params);
        
        // 出力ファイル設定
        SF_INFO sfinfo_out = sfinfo_in;
        sfinfo_out.samplerate = OUTPUT_SR;
        sfinfo_out.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT; // 32-bit float出力
        
        SNDFILE* outfile = sf_open(output_filename.c_str(), SFM_WRITE, &sfinfo_out);
        if (!outfile) {
            throw std::runtime_error("Could not create output file.");
        }
        
        // ダウンサンプラー
        Resampler downsampler(sfinfo_in.channels, (double)OUTPUT_SR / PROCESSING_SR);
        
        std::cout << "\nProcessing..." << std::endl;
        
        std::vector<float> input_buffer(BUFFER_SIZE * sfinfo_in.channels);
        sf_count_t total_frames_processed = 0;
        sf_count_t frames_read;
        
        // 処理前の音質分析
        std::vector<float> analysis_buffer;
        sf_count_t analysis_frames = std::min(static_cast<sf_count_t>(sfinfo_in.samplerate * 10), sfinfo_in.frames); // 最初の10秒
        analysis_buffer.resize(analysis_frames * sfinfo_in.channels);
        sf_seek(infile, 0, SEEK_SET);
        sf_readf_float(infile, analysis_buffer.data(), analysis_frames);
        auto before_analysis = AudioAnalyzer::analyze(analysis_buffer, sfinfo_in.channels, sfinfo_in.samplerate);
        
        // ファイルを最初に戻す
        sf_seek(infile, 0, SEEK_SET);
        
        // メイン処理ループ
        while ((frames_read = sf_readf_float(infile, input_buffer.data(), BUFFER_SIZE)) > 0) {
            input_buffer.resize(frames_read * sfinfo_in.channels);
            
            // アップサンプリング
            auto upsampled = upsampler.process(input_buffer, frames_read);
            
            // エフェクト処理
            processor.processBuffer(upsampled);
            
            // ダウンサンプリング
            auto final_output = downsampler.process(upsampled, upsampled.size() / sfinfo_in.channels);
            
            // 出力
            sf_writef_float(outfile, final_output.data(), final_output.size() / sfinfo_in.channels);
            
            total_frames_processed += frames_read;
            
            // 進捗表示
            if (total_frames_processed % (sfinfo_in.samplerate * 5) == 0) { // 5秒ごと
                double progress = (double)total_frames_processed / sfinfo_in.frames * 100.0;
                std::cout << "Progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::endl;
            }
        }
        
        sf_close(infile);
        sf_close(outfile);
        
        // 処理後の音質分析
        SNDFILE* analysis_file = sf_open(output_filename.c_str(), SFM_READ, &sfinfo_out);
        analysis_buffer.resize(analysis_frames * sfinfo_out.channels);
        sf_readf_float(analysis_file, analysis_buffer.data(), std::min(analysis_frames, sfinfo_out.frames));
        sf_close(analysis_file);
        
        auto after_analysis = AudioAnalyzer::analyze(analysis_buffer, sfinfo_out.channels, sfinfo_out.samplerate);
        
        // 結果表示
        std::cout << "\n=== Processing Complete ===" << std::endl;
        std::cout << "Output: " << output_filename << std::endl;
        std::cout << "\n=== Audio Analysis ===" << std::endl;
        std::cout << "Before -> After:" << std::endl;
        std::cout << "Peak Level: " << std::fixed << std::setprecision(2) 
                  << before_analysis.peak_db << "dB -> " << after_analysis.peak_db << "dB" << std::endl;
        std::cout << "RMS Level: " << before_analysis.rms_db << "dB -> " << after_analysis.rms_db << "dB" << std::endl;
        std::cout << "Dynamic Range: " << before_analysis.dynamic_range << "dB -> " << after_analysis.dynamic_range << "dB" << std::endl;
        
        if (sfinfo_in.channels == 2) {
            std::cout << "Stereo Width: " << std::fixed << std::setprecision(3)
                      << before_analysis.stereo_width << " -> " << after_analysis.stereo_width << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        sf_close(infile);
        return 1;
    }
    
    return 0;
}

// --- プリセット読み込み関数 ---
json loadPreset(const std::string& preset_name) {
    json params;
    
    if (preset_name == "master") {
        params = {
            {"spectral_gate", {{"enabled", false}}},
            {"multiband_comp", {
                {"enabled", true},
                {"bands", json::array({
                    {{"threshold_db", -15.0}, {"ratio", 2.5}, {"attack_ms", 10.0}, {"release_ms", 120.0}, {"makeup_gain_db", 1.5}, {"enabled", true}},
                    {{"threshold_db", -10.0}, {"ratio", 2.0}, {"attack_ms", 5.0}, {"release_ms", 80.0}, {"makeup_gain_db", 1.0}, {"enabled", true}},
                    {{"threshold_db", -8.0}, {"ratio", 1.8}, {"attack_ms", 3.0}, {"release_ms", 60.0}, {"makeup_gain_db", 1.5}, {"enabled", true}},
                    {{"threshold_db", -6.0}, {"ratio", 1.5}, {"attack_ms", 2.0}, {"release_ms", 40.0}, {"makeup_gain_db", 2.0}, {"enabled", true}}
                })}
            }},
            {"analog_saturation", {{"drive", 0.8}, {"mix", 0.15}, {"type", "tube"}}},
            {"harmonic_enhancer", {{"enabled", true}, {"drive", 0.4}, {"even_harmonics", 0.2}, {"odd_harmonics", 0.25}, {"mix", 0.2}}},
            {"stereo_enhancer", {{"enabled", true}, {"width", 1.15}, {"bass_mono_freq", 150.0}}},
            {"reverb", {{"enabled", false}}},
            {"final_limiter", {{"threshold_db", -0.2}, {"release_ms", 50.0}}}
        };
    } else if (preset_name == "vocal") {
        params = {
            {"spectral_gate", {{"enabled", true}, {"threshold_db", -45.0}}},
            {"harmonic_enhancer", {{"enabled", true}, {"drive", 0.6}, {"mix", 0.3}}},
            {"stereo_enhancer", {{"enabled", true}, {"width", 1.05}}},
            {"reverb", {{"enabled", true}, {"room_size", 0.2}, {"wet_level", 0.08}}}
        };
    } else if (preset_name == "instrumental") {
        params = {
            {"multiband_comp", {{"enabled", true}}},
            {"analog_saturation", {{"drive", 1.2}, {"mix", 0.25}, {"type", "tape"}}},
            {"stereo_enhancer", {{"enabled", true}, {"width", 1.3}}},
            {"reverb", {{"enabled", true}, {"room_size", 0.4}, {"wet_level", 0.12}}}
        };
    } else if (preset_name == "podcast") {
        params = {
            {"spectral_gate", {{"enabled", true}, {"threshold_db", -50.0}}},
            {"multiband_comp", {{"enabled", true}}},
            {"harmonic_enhancer", {{"enabled", true}, {"drive", 0.8}, {"mix", 0.4}}},
            {"stereo_enhancer", {{"enabled", false}}}
        };
    }
    
    // 外部設定ファイルがあれば上書き
    std::ifstream config_file("config_" + preset_name + ".json");
    if (config_file.is_open()) {
        json external_config;
        config_file >> external_config;
        params.update(external_config);
        std::cout << "Loaded external config: config_" << preset_name << ".json" << std::endl;
    }
    
    return params;
}
