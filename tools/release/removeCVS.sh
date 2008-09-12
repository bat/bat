#!/bin/sh

if [[ $1 == "" ]]; then
	EXCLUDE="CVS"
else
	EXCLUDE=$1
fi

for file in *; do
	if [ -d "$file" ]; then
		if [[ "$file" == "$EXCLUDE" ]]; then
			rm -rf $EXCLUDE
		else
			cd "$file"
			$0 $EXCLUDE
			cd ..
		fi
	fi
done
