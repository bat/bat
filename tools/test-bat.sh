#!/bin/bash -e

# Standard BAT build/test script (Linux)

./autogen.sh
./configure --with-cuba=download --enable-roostats --enable-parallel CXXFLAGS='-Wall -Wextra -pedantic -Werror'
make -j`nprocs` distcheck VERBOSE=1

echo "OK"
