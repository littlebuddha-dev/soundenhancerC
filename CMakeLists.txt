# ./CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(EnhancedSoundProcessor VERSION 2.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# コンパイラ最適化フラグ
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")
endif()

# 依存ライブラリを探す
find_package(PkgConfig REQUIRED)

# libsndfile
pkg_check_modules(SNDFILE REQUIRED sndfile)

# libsamplerate
pkg_check_modules(SAMPLERATE REQUIRED samplerate)

# FFTW3 (FFTベースのEQ用)
pkg_check_modules(FFTW3 fftw3f)

# nlohmann/json (ヘッダーオンリー)
find_path(NLOHMANN_JSON_INCLUDE_DIR nlohmann/json.hpp)

# 実行ファイル作成
add_executable(enhanced_sound_processor
    main.cpp
)

# インクルードディレクトリ
target_include_directories(enhanced_sound_processor PRIVATE
    ${SNDFILE_INCLUDE_DIRS}
    ${SAMPLERATE_INCLUDE_DIRS}
    ${NLOHMANN_JSON_INCLUDE_DIR}
)

if(FFTW3_FOUND)
    target_include_directories(enhanced_sound_processor PRIVATE ${FFTW3_INCLUDE_DIRS})
    target_compile_definitions(enhanced_sound_processor PRIVATE USE_FFTW3)
endif()

# ライブラリリンク
target_link_libraries(enhanced_sound_processor
    ${SNDFILE_LIBRARIES}
    ${SAMPLERATE_LIBRARIES}
)

if(FFTW3_FOUND)
    target_link_libraries(enhanced_sound_processor ${FFTW3_LIBRARIES})
endif()

# コンパイルフラグ
target_compile_options(enhanced_sound_processor PRIVATE
    ${SNDFILE_CFLAGS_OTHER}
    ${SAMPLERATE_CFLAGS_OTHER}
)

if(FFTW3_FOUND)
    target_compile_options(enhanced_sound_processor PRIVATE ${FFTW3_CFLAGS_OTHER})
endif()

# インストール設定
install(TARGETS enhanced_sound_processor DESTINATION bin)

# サンプル設定ファイルのインストール
install(FILES 
    config_master.json
    config_vocal.json
    config_instrumental.json
    config_podcast.json
    DESTINATION share/enhanced_sound_processor
)

---

# ./config_master.json
{
    "spectral_gate": {
        "enabled": false,
        "threshold_db": -60.0,
        "attack_ms": 5.0,
        "release_ms": 100.0
    },
    "multiband_comp": {
        "enabled": true,
        "bands": [
            {
                "name": "Bass",
                "freq_low": 20.0,
                "freq_high": 250.0,
                "threshold_db": -15.0,
                "ratio": 2.5,
                "attack_ms": 10.0,
                "release_ms": 120.0,
                "makeup_gain_db": 1.5,
                "enabled": true
            },
            {
                "name": "Low-Mid",
                "freq_low": 250.0,
                "freq_high": 2000.0,
                "threshold_db": -10.0,
                "ratio": 2.0,
                "attack_ms": 5.0,
                "release_ms": 80.0,
                "makeup_gain_db": 1.0,
                "enabled": true
            },
            {
                "name": "High-Mid",
                "freq_low": 2000.0,
                "freq_high": 8000.0,
                "threshold_db": -8.0,
                "ratio": 1.8,
                "attack_ms": 3.0,
                "release_ms": 60.0,
                "makeup_gain_db": 1.5,
                "enabled": true
            },
            {
                "name": "Treble",
                "freq_low": 8000.0,
                "freq_high": 20000.0,
                "threshold_db": -6.0,
                "ratio": 1.5,
                "attack_ms": 2.0,
                "release_ms": 40.0,
                "makeup_gain_db": 2.0,
                "enabled": true
            }
        ]
    },
    "analog_saturation": {
        "enabled": true,
        "drive": 0.8,
        "mix": 0.15,
        "type": "tube"
    },
    "harmonic_enhancer": {
        "enabled": true,
        "drive": 0.4,
        "even_harmonics": 0.2,
        "odd_harmonics": 0.25,
        "mix": 0.2,
        "freq_low": 100.0,
        "freq_high": 8000.0
    },
    "stereo_enhancer": {
        "enabled": true,
        "width": 1.15,
        "bass_mono_freq": 150.0
    },
    "reverb": {
        "enabled": false,
        "room_size": 0.3,
        "damping": 0.5,
        "wet_level": 0.1
    },
    "final_limiter": {
        "threshold_db": -0.2,
        "release_ms": 50.0
    }
}

---

# ./config_vocal.json
{
    "spectral_gate": {
        "enabled": true,
        "threshold_db": -45.0,
        "attack_ms": 3.0,
        "release_ms": 150.0
    },
    "multiband_comp": {
        "enabled": true,
        "bands": [
            {
                "freq_low": 20.0,
                "freq_high": 200.0,
                "threshold_db": -20.0,
                "ratio": 3.0,
                "attack_ms": 15.0,
                "release_ms": 100.0,
                "makeup_gain_db": 0.5,
                "enabled": true
            },
            {
                "freq_low": 200.0,
                "freq_high": 2000.0,
                "threshold_db": -12.0,
                "ratio": 2.2,
                "attack_ms": 5.0,
                "release_ms": 60.0,
                "makeup_gain_db": 2.0,
                "enabled": true
            },
            {
                "freq_low": 2000.0,
                "freq_high": 8000.0,
                "threshold_db": -8.0,
                "ratio": 1.8,
                "attack_ms": 2.0,
                "release_ms": 40.0,
                "makeup_gain_db": 3.0,
                "enabled": true
            },
            {
                "freq_low": 8000.0,
                "freq_high": 20000.0,
                "threshold_db": -6.0,
                "ratio": 1.5,
                "attack_ms": 1.0,
                "release_ms": 30.0,
                "makeup_gain_db": 2.5,
                "enabled": true
            }
        ]
    },
    "analog_saturation": {
        "enabled": true,
        "drive": 0.6,
        "mix": 0.1,
        "type": "tube"
    },
    "harmonic_enhancer": {
        "enabled": true,
        "drive": 0.6,
        "even_harmonics": 0.15,
        "odd_harmonics": 0.3,
        "mix": 0.3,
        "freq_low": 150.0,
        "freq_high": 12000.0
    },
    "stereo_enhancer": {
        "enabled": true,
        "width": 1.05,
        "bass_mono_freq": 200.0
    },
    "reverb": {
        "enabled": true,
        "room_size": 0.2,
        "damping": 0.6,
        "wet_level": 0.08
    },
    "final_limiter": {
        "threshold_db": -0.3,
        "release_ms": 30.0
    }
}

---

# ./build.sh
#!/bin/bash

# 依存関係のインストール (Ubuntu/Debian)
if command -v apt-get >/dev/null 2>&1; then
    echo "Installing dependencies for Ubuntu/Debian..."
    sudo apt-get update
    sudo apt-get install -y libsndfile1-dev libsamplerate0-dev nlohmann-json3-dev libfftw3-dev cmake build-essential pkg-config
fi

# 依存関係のインストール (macOS)
if command -v brew >/dev/null 2>&1; then
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

---

# ./README.md
# Enhanced Sound Processor

高品質なサウンドエンハンスメント処理を行うC++プログラムです。

## 特徴

### 高度な音響処理
- **マルチバンドコンプレッサー**: 4バンド構成で各帯域を独立制御
- **アナログサチュレーション**: チューブ、テープ、トランスフォーマーの特性を再現
- **ハーモニックエンハンサー**: 音楽的な倍音を付加して音に豊かさを追加
- **スペクトラルゲート**: 高度なノイズ除去
- **ステレオエンハンサー**: 立体感と空間の広がりを改善
- **高品質リバーブ**: 自然な残響効果

### 技術的特徴
- 192kHzでの内部処理による超高品質
- リンクウィッツ・ライリー4次クロスオーバーフィルター
- libsamplerateによる高品質リサンプリング
- 2パス処理による最適化されたワークフロー
- リアルタイム音質分析機能

## ビルド方法

```bash
# リポジトリをクローン
git clone <repository_url>
cd enhanced-sound-processor

# 依存関係のインストールとビルド
chmod +x build.sh
./build.sh
```

## 使用方法

### 基本的な使用法
```bash
# デフォルト（マスタリング）プリセット
./enhanced_sound_processor input.wav

# 特定のプリセットを指定
./enhanced_sound_processor input.wav vocal
./enhanced_sound_processor input.wav instrumental  
./enhanced_sound_processor input.wav podcast
```

### 利用可能なプリセット

#### master（デフォルト）
- 汎用的なマスタリング処理
- バランスの取れたダイナミクス処理
- 軽めのアナログ感とハーモニック強化

#### vocal
- ボーカル用に最適化
- ノイズゲート機能でクリーンな録音
- 中高域の明瞭度向上
- 自然なリバーブ効果

#### instrumental
- 楽器演奏用
- より強いアナログサチュレーション
- ワイドなステレオイメージ
- 豊かな空間表現

#### podcast
- 音声配信用
- 強力なノイズ除去
- 明瞭な音声強調
- ステレオ処理は控えめ

## 設定ファイル

各プリセットは対応するJSONファイルで詳細設定が可能です：

- `config_master.json`
- `config_vocal.json` 
- `config_instrumental.json`
- `config_podcast.json`

### 設定例
```json
{
    "harmonic_enhancer": {
        "enabled": true,
        "drive": 0.4,
        "even_harmonics": 0.2,
        "odd_harmonics": 0.25,
        "mix": 0.2
    },
    "stereo_enhancer": {
        "enabled": true,
        "width": 1.15,
        "bass_mono_freq": 150.0
    }
}
```

## 依存関係

- libsndfile (音声ファイルI/O)
- libsamplerate (高品質リサンプリング)
- nlohmann/json (設定ファイル処理)
- FFTW3 (オプション: FFTベースEQ用)
- CMake 3.15+
- C++17対応コンパイラ

## 対応フォーマット

- 入力: WAVファイル（各種ビット深度、サンプルレート）
- 出力: 32-bit float WAV @ 96kHz

## パフォーマンス

- 大容量バッファによる効率的処理
- マルチスレッド対応（将来の拡張予定）
- メモリ使用量の最適化
- リアルタイム処理能力

## 今後の改善予定

1. **FFTベースリニアフェーズEQ**: より正確な周波数特性制御
2. **AI支援音質分析**: 機械学習による最適パラメータ提案
3. **リアルタイム処理**: VST/AUプラグイン化
4. **高度なリバーブ**: コンボリューションリバーブ対応
5. **マルチスレッド処理**: 大容量ファイルの高速処理
