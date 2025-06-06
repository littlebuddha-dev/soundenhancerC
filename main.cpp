// ./main.cpp
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

// 3rd party libraries
#include <sndfile.h>
#include "json.hpp" // nlohmann/json

using json = nlohmann::json;

// --- Utilities ---
double db_to_linear(double db) {
    if (!std::isfinite(db)) return 1.0;
    return std::pow(10.0, db / 20.0);
}

float simple_clamp(float sample) {
    if (std::isnan(sample) || std::isinf(sample)) return 0.0f;
    return std::max(-1.0f, std::min(1.0f, sample));
}

// --- Filter Classes ---
class SimpleBiquad {
public:
    SimpleBiquad(std::string name = "Unnamed") : filter_name_(name) { reset(); }
    void reset() { a1 = 0.0; a2 = 0.0; b0 = 1.0; b1 = 0.0; b2 = 0.0; z1 = 0.0; z2 = 0.0; is_bypassed_ = false; }
    
    void set_lpf(double sr, double freq, double q) {
        q = std::max(0.1, q); freq = std::max(10.0, std::min(freq, sr / 2.2));
        double w0 = 2.0 * M_PI * freq / sr, cos_w0 = std::cos(w0), sin_w0 = std::sin(w0);
        double alpha = sin_w0 / (2.0 * q), a0 = 1.0 + alpha;
        b0 = (1.0 - cos_w0) / 2.0 / a0; b1 = (1.0 - cos_w0) / a0; b2 = b0;
        a1 = -2.0 * cos_w0 / a0; a2 = (1.0 - alpha) / a0;
    }
    void set_hpf(double sr, double freq, double q) {
        q = std::max(0.1, q); freq = std::max(10.0, std::min(freq, sr / 2.2));
        double w0 = 2.0 * M_PI * freq / sr, cos_w0 = std::cos(w0), sin_w0 = std::sin(w0);
        double alpha = sin_w0 / (2.0 * q), a0 = 1.0 + alpha;
        b0 = (1.0 + cos_w0) / 2.0 / a0; b1 = -(1.0 + cos_w0) / a0; b2 = b0;
        a1 = -2.0 * cos_w0 / a0; a2 = (1.0 - alpha) / a0;
    }
    void set_peaking(double sr, double freq, double q, double gain_db) {
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
        return simple_clamp(out);
    }
private:
    std::string filter_name_; bool is_bypassed_ = false;
    double a1, a2, b0, b1, b2, z1, z2;
};

// --- Processor Classes ---
class SimpleBBEProcessor {
public:
    void setup(double sr, const json& p) { if(!p.value("enabled",false))return; }
    float process(float in) { return in; }
};

class SimpleExciterProcessor {
public:
    SimpleExciterProcessor(): hpf_("Exciter_HPF"), interp_lpf_("InterpLPF"), decim_lpf_("DecimLPF") {}
    void setup(double sr, const json& p) {
        if(!p.value("enabled",false)) return; enabled_ = true;
        drive_ = p.value("drive", 5.0); mix_ = p.value("mix", 0.3); even_drive_ = p.value("even_drive", 0.5);
        hpf_.set_hpf(sr, p.value("crossover_freq", 6000.0), 0.707);
        interp_lpf_.set_lpf(sr * 2.0, sr / 2.2, 0.5); decim_lpf_.set_lpf(sr * 2.0, sr / 2.2, 0.5);
    }
    float process(float in) {
        if (!enabled_) return in;
        float high_freqs = hpf_.process(in);
        auto saturate = [&](float s) {
            float odd = std::tanh(s * drive_), even = (s * std::abs(s)) * even_drive_;
            return std::tanh(odd + even);
        };
        float s1 = saturate(interp_lpf_.process(high_freqs)); decim_lpf_.process(s1);
        float s2 = saturate(interp_lpf_.process(0.0f));
        float saturated_highs = decim_lpf_.process(s2);
        return simple_clamp(in * (1.0f - mix_) + saturated_highs * mix_);
    }
private:
    SimpleBiquad hpf_, interp_lpf_, decim_lpf_;
    float drive_, mix_, even_drive_; bool enabled_ = false;
};

// 改良版: 複数バンドのパラメトリックEQプロセッサー
class MultiParametricEQProcessor {
public:
    MultiParametricEQProcessor(): 
        peq1_("PEQ_Mid1"), peq2_("PEQ_Mid2"), peq3_("PEQ_Mid3"), 
        peq4_("PEQ_High1"), peq5_("PEQ_High2") {}
    
    void setup(double sr, const json& p) {
        if(!p.value("enabled",false)) return; 
        enabled_ = true;
        
        // 中域の補正用EQ（複数バンド）
        peq1_.set_peaking(sr, 2000.0, 1.2, 6.0);   // 2kHz +6dB (主要な中域補正)
        peq2_.set_peaking(sr, 4000.0, 1.0, 4.0);   // 4kHz +4dB (プレゼンス補正)
        peq3_.set_peaking(sr, 8000.0, 1.5, 3.0);   // 8kHz +3dB (高域中域補正)
        
        // 元の高域補正
        peq4_.set_peaking(sr, p.value("freq", 22000.0), p.value("q", 1.5), p.value("gain_db", 9.0));
        
        // 追加の超高域補正（必要に応じて）
        peq5_.set_peaking(sr, 16000.0, 2.0, 2.0);  // 16kHz +2dB
    }
    
    float process(float in) {
        if (!enabled_) return in;
        
        // 順次処理（周波数の低い順から）
        float out = peq1_.process(in);    // 2kHz補正
        out = peq2_.process(out);         // 4kHz補正
        out = peq3_.process(out);         // 8kHz補正
        out = peq5_.process(out);         // 16kHz補正
        out = peq4_.process(out);         // 25kHz補正
        
        return out;
    }
    
private:
    SimpleBiquad peq1_, peq2_, peq3_, peq4_, peq5_;
    bool enabled_ = false;
};

// --- Main Processing ---
int main(int argc, char *argv[]) {
    if (argc != 2) { std::cerr << "Usage: " << argv[0] << " <input_audio_file>" << std::endl; return 1; }
    
    std::string input_filename = argv[1];
    std::string output_filename = input_filename.substr(0, input_filename.rfind('.')) + "_final.wav";
    
    json params;
    std::ifstream f("params.json"); 
    if (f.is_open()) { 
        params = json::parse(f); 
        std::cout << "Parameters loaded from params.json" << std::endl;
    } else {
        std::cout << "Warning: params.json not found, using defaults" << std::endl;
    }
    
    SF_INFO sfinfo_in;
    SNDFILE* infile = sf_open(input_filename.c_str(), SFM_READ, &sfinfo_in);
    if (!infile) { std::cerr << "Error: Could not open input file." << std::endl; return 1; }
    
    SF_INFO sfinfo_out = sfinfo_in;
    sfinfo_out.samplerate *= 2; // Upsample output
    
    SNDFILE* outfile = sf_open(output_filename.c_str(), SFM_WRITE, &sfinfo_out);
    if (!outfile) { sf_close(infile); return 1; }
    
    std::vector<SimpleExciterProcessor> exciters(sfinfo_in.channels);
    std::vector<MultiParametricEQProcessor> peqs(sfinfo_in.channels);  // 改良版EQ使用
    std::vector<SimpleBiquad> resamplers(sfinfo_in.channels);

    for (int i = 0; i < sfinfo_in.channels; ++i) { 
        exciters[i].setup(sfinfo_in.samplerate, params["exciter"]); 
        peqs[i].setup(sfinfo_out.samplerate, params["peq"]); // PEQ runs at the high output sample rate
        resamplers[i].set_lpf(sfinfo_out.samplerate, sfinfo_in.samplerate / 2.2, 0.707);
    }
    
    const size_t BUF_SIZE = 4096;
    std::vector<float> in_buf(BUF_SIZE * sfinfo_in.channels), out_buf(BUF_SIZE * sfinfo_in.channels * 2);
    sf_count_t frames_read, total_frames = 0;
    
    std::cout << "Processing audio..." << std::endl;
    std::cout << "Input: " << sfinfo_in.samplerate << "Hz, " << sfinfo_in.channels << " channels" << std::endl;
    std::cout << "Output: " << sfinfo_out.samplerate << "Hz (upsampled)" << std::endl;
    
    while ((frames_read = sf_readf_float(infile, in_buf.data(), BUF_SIZE)) > 0) {
        // Step 1: Process with Exciter and upsample
        for (sf_count_t i = 0; i < frames_read; ++i) {
            for (int ch = 0; ch < sfinfo_in.channels; ++ch) {
                float sample = in_buf[i * sfinfo_in.channels + ch];
                float processed_sample = exciters[ch].process(sample);
                out_buf[(i * 2) * sfinfo_in.channels + ch] = processed_sample;
                out_buf[(i * 2 + 1) * sfinfo_in.channels + ch] = 0.0f;  // Zero-pad for upsampling
            }
        }
        
        // Step 2: Apply anti-aliasing filter and multi-band EQ
        for (sf_count_t i = 0; i < frames_read * 2; ++i) {
            for (int ch = 0; ch < sfinfo_in.channels; ++ch) {
                float sample = out_buf[i * sfinfo_in.channels + ch];
                sample = resamplers[ch].process(sample) * 2.0f;  // Anti-aliasing + gain compensation
                sample = peqs[ch].process(sample);  // Apply multi-band EQ after upsampling
                out_buf[i * sfinfo_in.channels + ch] = sample;
            }
        }
        
        sf_writef_float(outfile, out_buf.data(), frames_read * 2);
        total_frames += frames_read;
        
        // Progress indicator
        if (total_frames % (sfinfo_in.samplerate / 4) == 0) {
            std::cout << "." << std::flush;
        }
    }
    
    sf_close(infile); sf_close(outfile);
    std::cout << std::endl << "Processing complete. Output: " << output_filename << std::endl;
    std::cout << "Total frames processed: " << total_frames << std::endl;
    return 0;
}