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
    db = std::max(-60.0, std::min(60.0, db));
    return std::pow(10.0, db / 20.0);
}

// Simple clamp function
float simple_clamp(float sample) {
    if (std::isnan(sample) || std::isinf(sample)) return 0.0f;
    return std::max(-1.0f, std::min(1.0f, sample));
}

// --- Simplified Biquad Class ---
class SimpleBiquad {
public:
    SimpleBiquad(std::string name = "Unnamed") : filter_name_(name) {
        reset();
    }

    void reset() {
        a0 = 1.0; a1 = 0.0; a2 = 0.0;
        b0 = 1.0; b1 = 0.0; b2 = 0.0; // Passthrough
        z1 = 0.0; z2 = 0.0;
        is_bypassed_ = false;
    }

    void set_lpf(double sample_rate, double cutoff_freq, double q) {
        std::cout << "[" << filter_name_ << "] Setting LPF: " << cutoff_freq << "Hz" << std::endl;
        
        q = std::max(0.1, std::min(q, 10.0));
        cutoff_freq = std::max(10.0, std::min(cutoff_freq, sample_rate / 2.2));
        
        double w0 = 2.0 * M_PI * cutoff_freq / sample_rate;
        double sin_w0 = std::sin(w0);
        double cos_w0 = std::cos(w0);
        double alpha = sin_w0 / (2.0 * q);
        
        double temp_b1 = 1.0 - cos_w0;
        double temp_b0 = temp_b1 / 2.0;
        double temp_b2 = temp_b0;
        double temp_a0 = 1.0 + alpha;
        double temp_a1 = -2.0 * cos_w0;
        double temp_a2 = 1.0 - alpha;
        
        if (std::abs(temp_a0) < 1e-12) {
            std::cout << "[" << filter_name_ << "] Bypassing (a0 is too small)" << std::endl;
            is_bypassed_ = true;
            return;
        }
        
        b0 = temp_b0 / temp_a0;
        b1 = temp_b1 / temp_a0;
        b2 = temp_b2 / temp_a0;
        a1 = temp_a1 / temp_a0;
        a2 = temp_a2 / temp_a0;
        
        std::cout << "[" << filter_name_ << "] LPF setup complete" << std::endl;
    }

    void set_hpf(double sample_rate, double cutoff_freq, double q) {
        std::cout << "[" << filter_name_ << "] Setting HPF: " << cutoff_freq << "Hz" << std::endl;
        
        q = std::max(0.1, std::min(q, 10.0));
        cutoff_freq = std::max(10.0, std::min(cutoff_freq, sample_rate / 2.2));
        
        double w0 = 2.0 * M_PI * cutoff_freq / sample_rate;
        double sin_w0 = std::sin(w0);
        double cos_w0 = std::cos(w0);
        double alpha = sin_w0 / (2.0 * q);
        
        double temp_b0 = (1.0 + cos_w0) / 2.0;
        double temp_b1 = -(1.0 + cos_w0);
        double temp_b2 = temp_b0;
        double temp_a0 = 1.0 + alpha;
        double temp_a1 = -2.0 * cos_w0;
        double temp_a2 = 1.0 - alpha;
        
        if (std::abs(temp_a0) < 1e-12) {
            std::cout << "[" << filter_name_ << "] Bypassing (a0 is too small)" << std::endl;
            is_bypassed_ = true;
            return;
        }
        
        b0 = temp_b0 / temp_a0;
        b1 = temp_b1 / temp_a0;
        b2 = temp_b2 / temp_a0;
        a1 = temp_a1 / temp_a0;
        a2 = temp_a2 / temp_a0;
        
        std::cout << "[" << filter_name_ << "] HPF setup complete" << std::endl;
    }
    
    float process(float in_f) {
        if (is_bypassed_) {
            return simple_clamp(in_f);
        }
        
        if (std::isnan(in_f) || std::isinf(in_f)) {
            return 0.0f;
        }
        
        double in = static_cast<double>(in_f);
        
        // State check
        if (std::isnan(z1) || std::isinf(z1) || std::isnan(z2) || std::isinf(z2)) {
            z1 = z2 = 0.0;
            return simple_clamp(in_f);
        }

        // Filter calculation
        double out = b0 * in + z1;
        double new_z1 = b1 * in - a1 * out + z2;
        double new_z2 = b2 * in - a2 * out;
        
        // Result check
        if (std::isnan(out) || std::isinf(out) || std::isnan(new_z1) || 
            std::isinf(new_z1) || std::isnan(new_z2) || std::isinf(new_z2)) {
            z1 = z2 = 0.0;
            return simple_clamp(in_f);
        }

        z1 = new_z1;
        z2 = new_z2;
        
        return simple_clamp(static_cast<float>(out));
    }

private:
    std::string filter_name_;
    bool is_bypassed_ = false;
    double a0, a1, a2;
    double b0, b1, b2;
    double z1, z2;
};

class SimpleBBEProcessor {
public:
    SimpleBBEProcessor() : lpf_("BBE_LPF"), hpf_("BBE_HPF") {}

    void setup(double sample_rate, const json& params) {
        double crossover_low = params["crossover_freq_low"];
        double crossover_high = params["crossover_freq_high"];
        double low_boost_db = params["low_boost_db"];
        double high_boost_db = params["high_boost_db"];
        
        low_boost_ = db_to_linear(low_boost_db);
        high_boost_ = db_to_linear(high_boost_db);
        
        std::cout << "BBE: low=" << crossover_low << "Hz, high=" << crossover_high 
                  << "Hz, boost=" << low_boost_ << "/" << high_boost_ << std::endl;
        
        lpf_.set_lpf(sample_rate, crossover_low, 0.707);
        hpf_.set_hpf(sample_rate, crossover_high, 0.707);
    }

    float process(float in) {
        if (std::isnan(in) || std::isinf(in)) return 0.0f;
        
        float low_band = lpf_.process(in);
        float high_band = hpf_.process(in);
        float mid_band = in - low_band - high_band;
        
        if (std::isnan(low_band)) low_band = 0.0f;
        if (std::isnan(high_band)) high_band = 0.0f;
        if (std::isnan(mid_band)) mid_band = in;
        
        float result = (low_band * low_boost_) + mid_band + (high_band * high_boost_);
        return simple_clamp(result);
    }

private:
    SimpleBiquad lpf_, hpf_;
    float low_boost_ = 1.0f, high_boost_ = 1.0f;
};

class SimpleExciterProcessor {
public:
    SimpleExciterProcessor() : hpf_("Exciter_HPF") {}

    void setup(double sample_rate, const json& params) {
        double crossover_freq = params["crossover_freq"];
        drive_ = params["drive"];
        mix_ = params["mix"];
        
        std::cout << "Exciter: freq=" << crossover_freq << "Hz, drive=" << drive_ << ", mix=" << mix_ << std::endl;
        hpf_.set_hpf(sample_rate, crossover_freq, 0.707);
    }

    float process(float in) {
        if (std::isnan(in) || std::isinf(in)) return 0.0f;
        
        float high_freqs = hpf_.process(in);
        float drive_input = std::max(-10.0f, std::min(10.0f, high_freqs * drive_));
        if (std::isnan(drive_input)) drive_input = 0.0f;
        
        float saturated_highs = std::tanh(drive_input);
        if (std::isnan(saturated_highs)) saturated_highs = 0.0f;
        
        float result = in * (1.0f - mix_) + saturated_highs * mix_;
        return simple_clamp(result);
    }

private:
    SimpleBiquad hpf_;
    float drive_ = 0.5f;
    float mix_ = 0.3f;
};

// Function to detect silence
bool is_silence(const std::vector<float>& buffer, size_t frames, int channels) {
    for (size_t i = 0; i < frames * channels; ++i) {
        if (std::abs(buffer[i]) > 0.0001f) { // Treat even very small values as audio
            return false;
        }
    }
    return true;
}

// --- Main Processing ---
int main(int argc, char *argv[]) {
    if (argc != 2) { 
        std::cerr << "Usage: " << argv[0] << " <input_audio_file>" << std::endl; 
        return 1; 
    }
    
    std::string input_filename = argv[1];
    std::string output_filename;
    size_t dot_pos = input_filename.rfind('.');
    if (dot_pos != std::string::npos) { 
        output_filename = input_filename.substr(0, dot_pos) + "_processed" + input_filename.substr(dot_pos); 
    } else { 
        output_filename = input_filename + "_processed.wav"; 
    }
    
    // Default parameters
    json params = {
        {"bbe", {
            {"crossover_freq_low", 150.0},
            {"crossover_freq_high", 2500.0},
            {"low_band_delay_ms", 1.5},
            {"low_boost_db", 6.0},
            {"high_boost_db", 5.0}
        }},
        {"exciter", {
            {"crossover_freq", 4000.0},
            {"drive", 0.7},
            {"mix", 0.25}
        }}
    };
    
    try { 
        std::ifstream f("params.json"); 
        if (f.is_open()) {
            json file_params = json::parse(f);
            params = file_params;
            std::cout << "Loaded parameters from params.json." << std::endl;
        }
    }
    catch (...) { 
        std::cout << "Using default parameters." << std::endl;
    }
    
    // Open input file
    SF_INFO sfinfo_in;
    SNDFILE* infile = sf_open(input_filename.c_str(), SFM_READ, &sfinfo_in);
    if (!infile) { 
        std::cerr << "Error: Could not open input file" << std::endl; 
        return 1; 
    }
    
    std::cout << "=== File Information ===" << std::endl;
    std::cout << "Sample Rate: " << sfinfo_in.samplerate << " Hz" << std::endl;
    std::cout << "Channels: " << sfinfo_in.channels << std::endl;
    std::cout << "Frames: " << sfinfo_in.frames << std::endl;
    
    // Open output file
    SNDFILE* outfile = sf_open(output_filename.c_str(), SFM_WRITE, &sfinfo_in);
    if (!outfile) { 
        sf_close(infile); 
        return 1; 
    }
    
    // Initialize processors
    std::vector<SimpleBBEProcessor> bbe_processors(sfinfo_in.channels);
    std::vector<SimpleExciterProcessor> exciter_processors(sfinfo_in.channels);
    
    std::cout << "\n=== Initializing Processors ===" << std::endl;
    for (int i = 0; i < sfinfo_in.channels; ++i) { 
        std::cout << "Channel " << i << ":" << std::endl;
        bbe_processors[i].setup(sfinfo_in.samplerate, params["bbe"]); 
        exciter_processors[i].setup(sfinfo_in.samplerate, params["exciter"]); 
    }
    
    // Skip silence and find actual audio data
    const size_t SEARCH_BUFFER_SIZE = 4096;
    std::vector<float> search_buffer(SEARCH_BUFFER_SIZE * sfinfo_in.channels);
    sf_count_t frames_read;
    long long total_frames_searched = 0;
    bool found_audio = false;
    
    std::cout << "\n=== Searching for Audio Data ===" << std::endl;
    
    while ((frames_read = sf_readf_float(infile, search_buffer.data(), SEARCH_BUFFER_SIZE)) > 0) {
        if (!is_silence(search_buffer, frames_read, sfinfo_in.channels)) {
            found_audio = true;
            std::cout << "Audio data found! Position: " << total_frames_searched << " frames (" 
                      << (double)total_frames_searched / sfinfo_in.samplerate << " seconds)" << std::endl;
            
            // Display the first non-silent samples
            for (int i = 0; i < std::min((int)frames_read, 10); ++i) {
                for (int ch = 0; ch < sfinfo_in.channels; ++ch) {
                    float sample = search_buffer[i * sfinfo_in.channels + ch];
                    if (std::abs(sample) > 0.0001f) {
                        std::cout << "Frame " << (total_frames_searched + i) 
                                  << ", Ch " << ch << ": " << std::fixed << std::setprecision(6) << sample << std::endl;
                        break;
                    }
                }
            }
            break;
        }
        total_frames_searched += frames_read;
        
        if (total_frames_searched > 0 && total_frames_searched % (long long)(sfinfo_in.samplerate * 30) == 0) { // every 30 seconds
            std::cout << "Searching: " << (total_frames_searched / sfinfo_in.samplerate) << " seconds elapsed..." << std::endl;
        }
    }
    
    if (!found_audio) {
        std::cout << "Audio data not found. The entire file is silent." << std::endl;
        sf_close(infile);
        sf_close(outfile);
        return 1;
    }
    
    // Process the found audio data
    std::cout << "\n=== Starting Audio Processing ===" << std::endl;
    
    long long processed_frames = 0;
    int sample_count = 0;
    
    do {
        for (sf_count_t i = 0; i < frames_read; ++i) {
            for (int ch = 0; ch < sfinfo_in.channels; ++ch) {
                float sample = search_buffer[i * sfinfo_in.channels + ch];
                
                // Detailed display for the first few samples
                if (sample_count < 20 && std::abs(sample) > 0.0001f) {
                    std::cout << "\nSample #" << sample_count << " (Frame " << (total_frames_searched + processed_frames + i)
                              << ", Ch " << ch << "):" << std::endl;
                    std::cout << "  Input: " << std::fixed << std::setprecision(6) << sample << std::endl;
                    
                    float bbe_output = bbe_processors[ch].process(sample);
                    std::cout << "  After BBE: " << bbe_output << std::endl;
                    
                    float final_output = exciter_processors[ch].process(bbe_output);
                    std::cout << "  Final: " << final_output << std::endl;
                    
                    search_buffer[i * sfinfo_in.channels + ch] = final_output;
                    sample_count++;
                } else {
                    // Normal processing
                    sample = bbe_processors[ch].process(sample);
                    sample = exciter_processors[ch].process(sample);
                    search_buffer[i * sfinfo_in.channels + ch] = sample;
                }
            }
        }
        
        // Final safety check
        for (size_t idx = 0; idx < (size_t)frames_read * sfinfo_in.channels; ++idx) {
            if (std::isnan(search_buffer[idx]) || std::isinf(search_buffer[idx])) {
                search_buffer[idx] = 0.0f;
            }
        }
        
        sf_writef_float(outfile, search_buffer.data(), frames_read);
        processed_frames += frames_read;
        
        // Progress display
        if (processed_frames > 0 && processed_frames % (long long)(sfinfo_in.samplerate * 10) == 0) {
            std::cout << "Processed: " << (double)processed_frames / sfinfo_in.samplerate << " seconds" << std::endl;
        }
        
        // Read the next buffer
        frames_read = sf_readf_float(infile, search_buffer.data(), SEARCH_BUFFER_SIZE);
        
    } while (frames_read > 0);
    
    sf_close(infile);
    sf_close(outfile);
    
    // Display results
    std::ifstream output_check(output_filename, std::ios::binary | std::ios::ate);
    if (output_check.is_open()) {
        std::streamsize output_size = output_check.tellg();
        std::cout << "\n=== Results ===" << std::endl;
        std::cout << "Silent part: " << (double)total_frames_searched / sfinfo_in.samplerate << " seconds" << std::endl;
        std::cout << "Processed audio: " << (double)processed_frames / sfinfo_in.samplerate << " seconds" << std::endl;
        std::cout << "Output file size: " << (output_size / 1024.0 / 1024.0) << " MB" << std::endl;
        output_check.close();
    }
    
    std::cout << "Processing complete." << std::endl;
    return 0;
}