// ./main.cpp
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstdio> // For std::remove
#include <stdexcept>

// 3rd party libraries
#include <sndfile.h>
#include <samplerate.h> // libsamplerate header
#include "json.hpp" // nlohmann/json

using json = nlohmann::json;

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
    
    // Non-copyable
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


// --- Utilities ---
double db_to_linear(double db) {
    if (!std::isfinite(db)) return 1.0;
    return std::pow(10.0, db / 20.0);
}

// --- Filter Classes ---
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
        return out;
    }
private:
    std::string filter_name_; bool is_bypassed_ = false;
    double a1, a2, b0, b1, b2, z1, z2;
};

// --- Processor Classes ---
class SimpleExciterProcessor {
public:
    SimpleExciterProcessor(): hpf_("Exciter_HPF") {}
    void setup(double sr, const json& p) {
        if(!p.value("enabled",false)) return; enabled_ = true;
        drive_ = p.value("drive", 5.0); mix_ = p.value("mix", 0.3); even_drive_ = p.value("even_drive", 0.5);
        hpf_.set_hpf(sr, p.value("crossover_freq", 6000.0), 0.707);
    }
    float process(float in) {
        if (!enabled_) return in;
        float high_freqs = hpf_.process(in);
        auto saturate = [&](float s) {
            float odd = std::tanh(s * drive_);
            float even = (s * s - (1.0f/3.0f)) * even_drive_;
            return std::tanh(odd + even);
        };
        float saturated_highs = saturate(high_freqs);
        return in * (1.0f - mix_) + saturated_highs * mix_;
    }
private:
    SimpleBiquad hpf_;
    float drive_, mix_, even_drive_; bool enabled_ = false;
};


class MultiParametricEQProcessor {
public:
    // ★ 高域ブースト用のEQを撤廃
    MultiParametricEQProcessor(): 
        peq1_("PEQ_Mid1"), peq2_("PEQ_Mid2"), peq3_("PEQ_Mid3"), 
        lpf_final_("LPF_FinalCut") {}
    
    void setup(double sr, const json& p) {
        if(!p.value("enabled",false)) return; 
        enabled_ = true;
        
        // 中域のEQは維持
        peq1_.set_peaking(sr, 2000.0, 1.2, 3.0); 
        peq2_.set_peaking(sr, 4000.0, 1.0, 4.0);
        peq3_.set_peaking(sr, 8000.0, 1.5, 3.0);
        
        // ★ 高域をブーストするEQ (peq4, peq5) を削除

        // ★★★ LPFを調整して自然な減衰カーブを作る ★★★
        lpf_final_.set_lpf(sr, 35000.0, 0.707);
    }
    
    void process_inplace(std::vector<float>& buffer, int channels) {
        if (!enabled_) return;
        size_t frame_count = buffer.size() / channels;
        for (size_t i = 0; i < frame_count; ++i) {
            for (int ch = 0; ch < channels; ++ch) {
                float sample = buffer[i * channels + ch];
                sample = peq1_.process(sample);
                sample = peq2_.process(sample);
                sample = peq3_.process(sample);
                // ★ 高域EQの処理を削除
                sample = lpf_final_.process(sample);
                buffer[i * channels + ch] = sample;
            }
        }
    }
private:
    // ★ 高域ブースト用のEQを撤廃
    SimpleBiquad peq1_, peq2_, peq3_, lpf_final_;
    bool enabled_ = false;
};

// --- Main Processing ---
int main(int argc, char *argv[]) {
    if (argc != 2) { std::cerr << "Usage: " << argv[0] << " <input_audio_file>" << std::endl; return 1; }
    
    std::string input_filename = argv[1];
    std::string base_filename = input_filename.substr(0, input_filename.rfind('.'));
    std::string tmp_filename = base_filename + "_tmp.wav";
    std::string final_filename = base_filename + "_final.wav";
    
    json params;
    std::ifstream f("params.json");
    if (f.is_open()) {
        params = json::parse(f, nullptr, false);
        if (params.is_discarded()) {
             std::cerr << "Error parsing params.json. Using defaults." << std::endl;
             params = json({});
        } else {
            std::cout << "Parameters loaded from params.json." << std::endl;
        }
    } else {
        std::cout << "Warning: params.json not found, using defaults" << std::endl;
    }
    
    SF_INFO sfinfo_in;
    SNDFILE* infile = sf_open(input_filename.c_str(), SFM_READ, &sfinfo_in);
    if (!infile) {
        std::cerr << "Error: Could not open input file." << std::endl; return 1;
    }

    if ((sfinfo_in.format & SF_FORMAT_TYPEMASK) != SF_FORMAT_WAV) {
        std::cerr << "Error: This program only supports WAV files." << std::endl;
        sf_close(infile); return 1;
    }

    const int TMP_SR = 192000;
    const int FINAL_SR = 96000;
    const size_t BUF_SIZE = 4096;

    try {
        // =========================================================================
        //  第1パス：アップサンプル、エフェクト処理、一時ファイルへの書き出し
        // =========================================================================
        std::cout << "\n--- Pass 1: Upsampling & Applying effects to 32-bit float temp file ---" << std::endl;
        
        SF_INFO sfinfo_tmp = sfinfo_in;
        sfinfo_tmp.samplerate = TMP_SR;
        sfinfo_tmp.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

        SNDFILE* tmp_outfile = sf_open(tmp_filename.c_str(), SFM_WRITE, &sfinfo_tmp);
        if (!tmp_outfile) {
            throw std::runtime_error("Could not create temp file.");
        }

        Resampler upsampler(sfinfo_in.channels, (double)TMP_SR / sfinfo_in.samplerate);
        std::vector<SimpleExciterProcessor> exciters(sfinfo_in.channels);
        MultiParametricEQProcessor eq_processor; // 1 instance for all channels
        
        for (int i = 0; i < sfinfo_in.channels; ++i) {
            exciters[i].setup(TMP_SR, params.value("exciter", json({})));
        }
        eq_processor.setup(TMP_SR, params.value("peq", json({})));
        
        std::vector<float> in_buf(BUF_SIZE * sfinfo_in.channels);
        sf_count_t frames_read;
        
        while ((frames_read = sf_readf_float(infile, in_buf.data(), BUF_SIZE)) > 0) {
            in_buf.resize(frames_read * sfinfo_in.channels);
            std::vector<float> upsampled_buf = upsampler.process(in_buf, frames_read);

            long upsampled_frames = upsampled_buf.size() / sfinfo_in.channels;
            for(long i = 0; i < upsampled_frames; ++i) {
                for(int ch = 0; ch < sfinfo_in.channels; ++ch) {
                    upsampled_buf[i * sfinfo_in.channels + ch] = exciters[ch].process(upsampled_buf[i * sfinfo_in.channels + ch]);
                }
            }
            eq_processor.process_inplace(upsampled_buf, sfinfo_in.channels);

            sf_writef_float(tmp_outfile, upsampled_buf.data(), upsampled_frames);
        }
        sf_close(infile);
        sf_close(tmp_outfile);

        // =========================================================================
        //  第2パス：ノーマライズとダウンサンプル、最終ファイルへの書き出し
        // =========================================================================
        std::cout << "\n--- Pass 2: Normalizing and downsampling to final format ---" << std::endl;

        SNDFILE* tmp_infile = sf_open(tmp_filename.c_str(), SFM_READ, &sfinfo_tmp);
        if (!tmp_infile) {
            throw std::runtime_error("Could not open temp file for reading.");
        }
        
        float max_peak = 0.0f;
        std::vector<float> read_buf(BUF_SIZE * sfinfo_tmp.channels);
        sf_count_t read_count;
        while ((read_count = sf_readf_float(tmp_infile, read_buf.data(), BUF_SIZE)) > 0) {
            for (sf_count_t i = 0; i < read_count * sfinfo_tmp.channels; ++i) {
                if (std::abs(read_buf[i]) > max_peak) {
                    max_peak = std::abs(read_buf[i]);
                }
            }
        }
        std::cout << "Peak value found: " << max_peak << std::endl;

        float norm_gain = (max_peak > 0.0f) ? (float)(std::pow(10.0, -0.1 / 20.0) / max_peak) : 1.0f;
        std::cout << "Normalization gain: " << norm_gain << std::endl;

        sf_seek(tmp_infile, 0, SEEK_SET);

        SF_INFO sfinfo_final = sfinfo_tmp;
        sfinfo_final.samplerate = FINAL_SR;
        sfinfo_final.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;

        SNDFILE* final_outfile = sf_open(final_filename.c_str(), SFM_WRITE, &sfinfo_final);
        if (!final_outfile) {
            throw std::runtime_error("Could not create final output file.");
        }
        
        Resampler downsampler(sfinfo_final.channels, (double)FINAL_SR / TMP_SR);
        
        while ((read_count = sf_readf_float(tmp_infile, read_buf.data(), BUF_SIZE)) > 0) {
            read_buf.resize(read_count * sfinfo_tmp.channels);
            for (sf_count_t i = 0; i < read_count * sfinfo_tmp.channels; ++i) {
                read_buf[i] *= norm_gain;
            }
            std::vector<float> final_buf = downsampler.process(read_buf, read_count);
            sf_writef_float(final_outfile, final_buf.data(), final_buf.size() / sfinfo_final.channels);
        }

        sf_close(tmp_infile);
        sf_close(final_outfile);

        std::remove(tmp_filename.c_str());

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        sf_close(infile); // Close infile if it's still open
        // Clean up temporary file if it exists
        std::ifstream temp_file(tmp_filename);
        if (temp_file.good()) {
            temp_file.close();
            std::remove(tmp_filename.c_str());
        }
        return 1;
    }

    std::cout << "\nProcessing complete. Output: " << final_filename << std::endl;
    return 0;
}
