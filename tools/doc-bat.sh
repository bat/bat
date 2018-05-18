#!/bin/bash -e

# Standard BAT build and deploy docs script (Linux)

./autogen.sh
# don't enable parallel to test that serial code builds as well
./configure --with-cuba=download --enable-roostats --enable-debug CXXFLAGS='-Werror'
# build html reference guide, manual in html and pdf
make -j`nproc` && make -C doc/manual/ pdf

echo "OK"
