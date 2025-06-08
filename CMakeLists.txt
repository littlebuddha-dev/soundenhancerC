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
    params.json
    config_master.json
    config_vocal.json
    config_instrumental.json
    config_podcast.json
    DESTINATION share/${PROJECT_NAME}
)