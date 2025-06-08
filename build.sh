#!/bin/bash

# 依存関係のインストール (Ubuntu/Debian)
if command -v apt-get >/dev/null 2>&1;
then
    echo "Installing dependencies for Ubuntu/Debian..."
    sudo apt-get update
    sudo apt-get install -y libsndfile1-dev libsamplerate0-dev nlohmann-json3-dev libfftw3-dev cmake build-essential pkg-config
fi

# 依存関係のインストール (macOS)
if command -v brew >/dev/null 2>&1;
then
    echo "Installing dependencies for macOS..."
    brew install libsndfile libsamplerate nlohmann-json fftw cmake pkg-config
fi

# ビルドディレクトリ作成
mkdir -p build
cd build

# CMake設定
cmake .. -DCMAKE_BUILD_TYPE=Release

# ビルド実行
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo "Build complete! Executable: ./enhanced_sound_processor"
echo ""
echo "Usage examples:"
echo "  ./enhanced_sound_processor input.wav"
echo "  ./enhanced_sound_processor input.wav master"
echo "  ./enhanced_sound_processor input.wav vocal"
echo "  ./enhanced_sound_processor input.wav instrumental"
echo "  ./enhanced_sound_processor input.wav podcast"