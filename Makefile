#
# magnon dynamics module for Takin
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2022
# @license GPLv2 (see 'LICENSE' file)
#

mingw_build = 0
debug_build = 0
strip_bins = 1


# -----------------------------------------------------------------------------
# setup
# -----------------------------------------------------------------------------
ifneq ($(mingw_build), 1)
	ifeq ("$(CXX)", "")
		CXX = g++
	endif

	SYSINCS = -I/usr/local/include \
		-I/usr/include/lapacke -I/usr/local/opt/lapack/include \
		-I/usr/include/qt5 -I/usr/include/x86_64-linux-gnu/qt5/ \
		-I/usr/local/include/Minuit2 \
		#-I/usr/local/Cellar/qt/5.15.0/include \
		#-I/home/tw/build/boost_1_73_0
	LIBDIRS = -L/usr/local/opt/lapack/lib -L/usr/local/lib
else
	CXX = x86_64-w64-mingw32-g++

	SYSINCS = -I/usr/x86_64-w64-mingw32/sys-root/mingw/include \
		-I/usr/x86_64-w64-mingw32/sys-root/mingw/include/qt5 \
		-I/usr/x86_64-w64-mingw32/sys-root/mingw/include/Minuit2
	LIBDIRS = -L/usr/x86_64-w64-mingw32/sys-root/mingw/bin/
endif


ifneq ($(debug_build), 1)
	OPT = -O2 -march=native

	ifeq ($(strip_bins), 1)
		STRIP = strip
	else
		STRIP = echo -e "Not stripping"
	endif
else
	OPT = -g -ggdb -Wall -Wextra
	STRIP = echo -e "Not stripping"
endif


STD = -std=c++20
LIBDEFS = -fPIC -DUSE_LAPACK=1
INCS = -Isrc -Iext -Iext/takin $(SYSINCS)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# meta rules
# -----------------------------------------------------------------------------
.PHONY: all clean

all: prepare \
	lib/magnonmod.so

clean:
	find . -name "*.o" -exec rm -fv {} \;
	rm -rfv bin/
	rm -rfv lib/

prepare:
	mkdir -p bin/
	mkdir -p lib/
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Takin plugin modules
# -----------------------------------------------------------------------------
lib/magnonmod.so: src/magnonmod.o \
		ext/takin/tools/monteconvo/sqwbase.o \
		ext/tlibs2/libs/log.o ext/tlibs/log/log.o
	@echo "Linking Takin module $@..."
	$(CXX) $(STD) $(OPT) $(LIBDIRS) $(LIBDEFS) -shared -o $@ $+ \
		-llapacke
	$(STRIP) $@
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# general rules
# -----------------------------------------------------------------------------
%.o: %.cpp
	@echo "Compiling $< -> $@..."
	$(CXX) $(STD) $(OPT) $(INCS) $(LIBDEFS) -c $< -o $@
# -----------------------------------------------------------------------------
