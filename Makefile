# ./Makefile
# 現在のmain.cppファイル用の簡易Makefile

CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra
TARGET = sound_processor
SOURCE = main.cpp

# 依存ライブラリ
LIBS = -lsndfile -lsamplerate

# pkg-configを使ってフラグを自動取得
CXXFLAGS += $(shell pkg-config --cflags sndfile)
CXXFLAGS += $(shell pkg-config --cflags samplerate)
LIBS += $(shell pkg-config --libs sndfile)
LIBS += $(shell pkg-config --libs samplerate)

# nlohmann/jsonの場所を探す（システムにインストールされている場合）
JSON_INCLUDE = $(shell find /usr/include /usr/local/include -name "nlohmann" -type d 2>/dev/null | head -1)
ifneq ($(JSON_INCLUDE),)
    CXXFLAGS += -I$(dir $(JSON_INCLUDE))
endif

# デバッグビルド用
DEBUG_FLAGS = -g -DDEBUG -O0
RELEASE_FLAGS = -DNDEBUG -O3 -march=native

.PHONY: all clean debug release install

all: release

release:
	$(CXX) $(CXXFLAGS) $(RELEASE_FLAGS) $(SOURCE) -o $(TARGET) $(LIBS)
	@echo "Release build complete: ./$(TARGET)"

debug:
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) $(SOURCE) -o $(TARGET)_debug $(LIBS)
	@echo "Debug build complete: ./$(TARGET)_debug"

clean:
	rm -f $(TARGET) $(TARGET)_debug

install: release
	sudo cp $(TARGET) /usr/local/bin/
	sudo cp params.json /usr/local/share/$(TARGET)/
	@echo "Installed to /usr/local/bin/$(TARGET)"

# 依存関係チェック
check-deps:
	@echo "Checking dependencies..."
	@pkg-config --exists sndfile || (echo "ERROR: libsndfile not found. Install with: sudo apt-get install libsndfile1-dev" && exit 1)
	@pkg-config --exists samplerate || (echo "ERROR: libsamplerate not found. Install with: sudo apt-get install libsamplerate0-dev" && exit 1)
	@test -f /usr/include/nlohmann/json.hpp || test -f /usr/local/include/nlohmann/json.hpp || (echo "ERROR: nlohmann/json not found. Install with: sudo apt-get install nlohmann-json3-dev" && exit 1)
	@echo "All dependencies found!"

# ヘルプ
help:
	@echo "Available targets:"
	@echo "  all          - Build release version (default)"
	@echo "  release      - Build optimized release version"
	@echo "  debug        - Build debug version"
	@echo "  clean        - Remove built files"
	@echo "  install      - Install to system (requires sudo)"
	@echo "  check-deps   - Check if dependencies are installed"
	@echo "  help         - Show this help"
	@echo ""
	@echo "Usage examples:"
	@echo "  make check-deps"
	@echo "  make release"
	@echo "  ./$(TARGET) input.wav"

# 依存ライブラリのインストールヘルプ
install-deps-ubuntu:
	sudo apt-get update
	sudo apt-get install -y libsndfile1-dev libsamplerate0-dev nlohmann-json3-dev build-essential pkg-config

install-deps-macos:
	brew install libsndfile libsamplerate nlohmann-json pkg-config