SCRIPT_SOURCE = blitzDrt.c++
SCRIPT_NAME = productiveDrt
HEADER_NAME = sourcecode.h

PREFIX = /usr/local

CMAKE_DIR = cmake-3.19.6-Linux-x86_64
CMAKE_INSTALLER = https://github.com/Kitware/CMake/releases/download/v3.19.6/$(CMAKE_DIR).sh

RM = rm -vrf

PATH := $(PWD)/$(CMAKE_DIR)/bin:$(PATH)
export PATH

# BEGIN-EVAL makefile-parser --make-help Makefile

help:
	@echo ""
	@echo "  Targets"
	@echo ""
	@echo "    deps-ubuntu    Install Ubuntu 18.04 packages"
	@echo "    check          Check whether dependencies are available"
	@echo "    install-cmake  Install $(notdir $(CMAKE_INSTALLER))"
	@echo "    install-blitz  build and install blitz"
	@echo "    build          build the binary"
	@echo ""
	@echo "  Variables"
	@echo ""

# END-EVAL

# Install Ubuntu 18.04 packages
deps-ubuntu:
	sudo apt-get install libboost-program-options-dev libfftw3-dev

# Check whether dependencies are available
check:
	@pkg-config fftw3 || echo "fftw3 missing, install with 'make deps-ubuntu'"
	@if ! pkg-config blitz;then \
		echo "blitz not installed. Install with 'make install-blitz'"; \
		cmake_ver=`cmake --version |grep -o '[0-9].*'`; \
		case "$$cmake_ver" in \
			1*|2*|3.1.*|3.2.*|3.3.*|3.4.*|3.5.*|3.6.*|3.7.*|3.8.*|3.9.*|3.10.*|3.11.*) \
				echo "cmake $$cmake_ver is too old, need 3.12.0+ to build blitz. Try make install-cmake" ;; \
		esac; \
	fi

# Install $(notdir $(CMAKE_INSTALLER)). Needed to build blitz if local cmake is too old
install-cmake: $(CMAKE_DIR).sh
	sh ./$(CMAKE_DIR).sh

$(CMAKE_DIR).sh:
	wget $(CMAKE_INSTALLER)

# build and install blitz
install-blitz:
	git submodule update
	cd repo/blitz && \
	mkdir -p build && \
	cd build && \
	cmake .. && \
	make lib && \
	sudo make install

# build the binary
build: $(SCRIPT_NAME)

install: build
	sudo install $(SCRIPT_NAME) $(PREFIX)/bin

uninstall:
	sudo $(RM) $(PREFIX)/bin/$(SCRIPT_NAME)

$(SCRIPT_NAME): $(HEADER_NAME)
	g++ -o $@ \
		-O3 \
		-march=native -mtune=native \
		$(SCRIPT_SOURCE) \
		-Wall -Wextra \
		-lfftw3 \
		-lboost_program_options \
		`pkg-config --cflags --libs blitz` \
		`pkg-config --cflags --libs fftw3` \
		`pkg-config --cflags --libs Magick++`
	strip $@

$(HEADER_NAME):
	./c2macro.pl $(SCRIPT_SOURCE) > $@

# remove built files
clean:
	$(RM) $(SCRIPT_NAME)
	$(RM) $(HEADER_NAME)
	$(RM) $(CMAKE_DIR)
	$(RM) $(CMAKE_DIR).sh
	cd repo/blitz; git checkout .; git clean -df;
