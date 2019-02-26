#!/usr/bin/env bash

export THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export PROJECT_ROOT=$THIS_DIR/..

export CC=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-gcc
export CXX=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-g++

# debug flags:
export CFLAGS="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -I${CONDA_PREFIX}/include"
export DEBUG_CXXFLAGS="-fvisibility-inlines-hidden -std=c++17 -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -I${CONDA_PREFIX}/include"


mkdir -p $THIS_DIR/build
mkdir -p $THIS_DIR/lib
cd $THIS_DIR/build

cmake $PROJECT_ROOT -DCMAKE_INSTALL_PREFIX=$THIS_DIR -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_CXX_FLAGS="${DEBUG_CXXFLAGS}" -DCMAKE_C_FLAGS="${DEBUG_CFLAGS}"
make VERBOSE=1
make install
mv *.so ../lib
