#!/bin/bash -e

# Standard BAT build/test script (Linux)

./autogen.sh
./configure --with-cuba=download --enable-roostats --enable-parallel
make -j`nprocs` distcheck VERBOSE=1 CXXFLAGS='-Wall -Wextra -pedantic -Werror'

echo "OK"
