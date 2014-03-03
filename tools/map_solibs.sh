#!/bin/bash

CLASS=$1
USERLIB=""
OBJECT=""
if [[ x"$3" == x ]]; then
	USERLIB=$2
	OBJECT=$USERLIB.so
else
	OBJECT=$2
	USERLIB=$3
fi

#from Axel: http://root.cern.ch/phpBB2/viewtopic.php?t=4778
SYMBOLSROOT=$(symbols=`nm -C $OBJECT | grep -E ' U T[[:alnum:]_]+::' | sed 's,^.* U \([[:alnum:]_]\+\)::.*$,(\1),'| sort | uniq | tr '\n' '|'`; (for so in $ROOTSYS/lib/*.so; do nm -C -D --defined-only $so | grep -E 'T ('$symbols')::' > /dev/null && echo $so; done )| sort | uniq )
SYMBOLSROOSTATS=$(symbols=`nm -C $OBJECT | grep -E ' U Roo[[:alnum:]_]+::' | sed 's,^.* U \([[:alnum:]_]\+\)::.*$,(\1),'| sort | uniq | tr '\n' '|'`; (for so in $ROOTSYS/lib/*.so; do nm -C -D --defined-only $so | grep -E 'T ('$symbols')::' > /dev/null && echo $so; done )| sort | uniq )
SYMBOLSBAT=$(symbols=`nm -C $OBJECT | grep -E ' U BC[[:alnum:]_]+::' | sed 's,^.* U \([[:alnum:]_]\+\)::.*$,(\1),'| sort | uniq | tr '\n' '|'`; (for so in $BATINSTALLDIR/lib/*.so; do nm -C -D --defined-only $so | grep -E 'T ('$symbols')::' > /dev/null && echo $so; done )| sort | uniq )

#remove $ROOTSYS
SSYMBOLSROOT=$(echo $SYMBOLSROOT | sed s="$ROOTSYS/lib/"==g)
SSYMBOLSROOSTATS=$(echo $SYMBOLSROOSTATS | sed s="$ROOTSYS/lib/"==g)
SSYMBOLSBAT=$(echo $SYMBOLSBAT | sed s="$BATINSTALLDIR/lib/"==g)

#now output rootmap
echo -n "Library.$CLASS:       $USERLIB.so "
for BATLIB in $SSYMBOLSBAT ; do
	if [[ "$BATLIB" == "$USERLIB.so" ]]; then
		true
	else
		echo -n "$BATLIB "
	fi
done
for ROOTLIB in $SSYMBOLSROOT ; do
    echo -n "$ROOTLIB "
done
for ROOSTATSLIB in $SSYMBOLSROOSTATS ; do
    echo -n "$ROOSTATSLIB "
done
echo ""
