This program is written in C++ to multiply BBE and exciter to wav files.

$ brew install libsndfile
$  brew install nlohmann-json

Download the json.hpp file from the link below and put it in the same hierarchy as main.cpp
https://github.com/nlohmann/json/blob/develop/single_include/nlohmann/json.hpp

# Install libsamplerate
$ brew install libsamplerate

# Create and move to build directory
$ mkdir build
$ cd build

# Run CMake to generate build files
$ cmake ..

# Execute compilation
$ make



$ apply_effects your_audiofile.wav
