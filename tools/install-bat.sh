#!/bin/bash -e

INSTALL_PREFIX="$1"

if [ -z "${INSTALL_PREFIX}" ] ; then
    echo "Syntax: $0 INSTALL_PREFIX" >&2
fi


SRCDIR="$PWD"
test -f configure || ./autogen.sh

BUILD_AREA="`mktemp -d -t build-bat-XXXXXXXX`"

echo "INFO: Build area: \"${BUILD_AREA}\"" >&2

(
    cd "${BUILD_AREA}"

    "${SRCDIR}/configure" \
        --prefix="$INSTALL_PREFIX" \
        --with-cuba=download \
        --enable-roostats \
        --enable-parallel

    make -j`nprocs` install
)

rm -rf "${BUILD_AREA}"

echo "Installation successful"
