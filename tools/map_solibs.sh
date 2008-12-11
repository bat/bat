#!/bin/bash

CLASS=$1
USERLIB=$2


#from Axel: http://root.cern.ch/phpBB2/viewtopic.php?t=4778
SYMBOLS=$(symbols=`nm -C $USERLIB.so | grep -E ' U T[[:alnum:]_]+::' | sed 's,^.* U \([[:alnum:]_]\+\)::.*$,(\1),'| sort | uniq | tr '\n' '|'`; (for so in $ROOTSYS/lib/root/*.so; do nm -C -D --defined-only $so | grep -E 'T ('$symbols')::' > /dev/null && echo $so; done )| sort | uniq )

#remove $ROOTSYS
SSYMBOLS=$(echo $SYMBOLS | sed s="$ROOTSYS/lib/"==g)

#now output rootmap
echo -n "Library.$CLASS:       $USERLIB.so "
for ROOTLIB in $SSYMBOLS ; do
    echo -n "$ROOTLIB "
done
echo ""
