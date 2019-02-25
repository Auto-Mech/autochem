#!/usr/bin/env bash
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_CXX_FLAGS="${DEBUG_CXXFLAGS}" -DCMAKE_C_FLAGS="${DEBUG_CFLAGS}"
make VERBOSE=1
make install
cd ..
mv build/*.so .
$PYTHON setup.py install
