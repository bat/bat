#!/bin/bash -e

CUBA_VERSION_FILE="cuba-version.txt"
LOCAL_CUBA_DIR="cuba"

if [ -f "${CUBA_VERSION_FILE}" ] ; then
	CUBA_VERSION=`cat "${CUBA_VERSION_FILE}" | head -n1`
else
	echo "Error: No file \"${CUBA_VERSION_FILE}\" here - must run in BAT top source directory." >&2
	exit 1
fi

CUBA_URL="http://www.feynarts.de/cuba/Cuba-${CUBA_VERSION}.tar.gz"


if [ -d "${LOCAL_CUBA_DIR}" ] ; then
	echo "Error: Directory ${LOCAL_CUBA_DIR} exists, please remove it first." >&2
	exit 1
fi


mkdir "${LOCAL_CUBA_DIR}"

if (which curl > /dev/null) ; then
	DL_COMMAND="curl -L"
elif (which wget > /dev/null) ; then
	DL_COMMAND="wget -O-"
else
	echo "Error: Need either \"curl\" or \"wget\", but neither seems to be on your PATH." >&2
	exit 1	
fi

${DL_COMMAND} "${CUBA_URL}" | tar --strip-components 1 -C cuba --strip=1 -x -z


(cd "${LOCAL_CUBA_DIR}" && CFLAGS="-fPIC -O3 -fomit-frame-pointer -ffast-math -Wall" ./configure && make)
