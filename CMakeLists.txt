cmake_minimum_required(VERSION 3.10)
project(AudioEffects)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- 外部ライブラリの検索 ---

# libsndfileライブラリを検索します。
# これにより、SNDFILE_INCLUDE_DIR と SNDFILE_LIBRARIES という変数が設定されます。
find_package(SndFile REQUIRED)

# nlohmann/jsonはヘッダーオンリーなので、ここでは特別な検索は不要です。


# --- 実行ファイルの定義 ---
add_executable(apply_effects main.cpp)


# --- 依存関係の設定 ---

# 1. libsndfileのヘッダーファイルへのパスを設定します。
#    これにより、コンパイラが <sndfile.h> を見つけられるようになります。
target_include_directories(apply_effects PRIVATE ${SNDFILE_INCLUDE_DIR})

# 2. カレントのソースディレクトリをインクルードパスに追加します。
#    これにより、コンパイラが "json.hpp" を見つけられるようになります。
#    ※必ず main.cpp と同じディレクトリに json.hpp を置いてください。
target_include_directories(apply_effects PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# 3. libsndfileライブラリ本体をリンクします。
#    エラーの原因となっていた SndFile::SndFile の代わりに、
#    より多くの環境で動作する変数 ${SNDFILE_LIBRARIES} を使用します。
target_link_libraries(apply_effects PRIVATE ${SNDFILE_LIBRARIES})