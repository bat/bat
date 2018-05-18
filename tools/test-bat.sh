#!/bin/bash -e

# Standard BAT build/test script (Linux)

./autogen.sh
./configure --with-cuba=download --enable-roostats --enable-parallel CXXFLAGS='-Wall -Wextra -pedantic -Werror'
# build html reference guide, manual in html and pdf
make -j`nproc` distcheck VERBOSE=1

echo "OK"
