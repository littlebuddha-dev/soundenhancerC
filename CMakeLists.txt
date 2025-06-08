# ./CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(sound_enhancer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- ライブラリの探索 (pkg-configに統一) ---
# find_package(Sndfile) の代わりに pkg-config を使用します。
find_package(PkgConfig REQUIRED)

# pkg-configを使ってsndfileライブラリを探索
pkg_check_modules(SNDFILE REQUIRED sndfile)

# pkg-configを使ってsamplerateライブラリを探索
pkg_check_modules(SAMPLERATE REQUIRED samplerate)

# --- 実行ファイルの定義 ---
add_executable(apply_effects main.cpp)

# --- インクルードディレクトリの追加 ---
# pkg-configが見つけたパスを使用
target_include_directories(apply_effects PRIVATE 
    ${SNDFILE_INCLUDE_DIRS}
    ${SAMPLERATE_INCLUDE_DIRS}
)

# --- ライブラリのリンク ---
# pkg-configが見つけたライブラリをリンク
target_link_libraries(apply_effects PRIVATE 
    ${SNDFILE_LIBRARIES}
    ${SAMPLERATE_LIBRARIES}
)
