#!/bin/bash

# Precommit hook for git that checks the code style.

# Check the wiki page at https://github.com/bat/bat/wiki/BAT-Code-Style,
# then activate the hook:
# cd .git/hooks/
# ln -s ../../tools/pre-commit-hook.sh pre-commit

OPTIONS="--options=tools/astylerc"

ASTYLE=$(which astyle)
if [ $? -ne 0 ]; then
	echo "[commit rejected] astyle not installed. Unable to check source file format policy." >&2
	exit 1
fi

# compare version strings a and b. return
#  9) a < b
# 10) a == b
# 11) a > b
# taken from http://stackoverflow.com/a/4024038/987623
do_version_check() {

    [ "$1" == "$2" ] && return 10

    ver1front=`echo $1 | cut -d "." -f -1`
    ver1back=`echo $1 | cut -d "." -f 2-`

    ver2front=`echo $2 | cut -d "." -f -1`
    ver2back=`echo $2 | cut -d "." -f 2-`

    if [ "$ver1front" != "$1" ] || [ "$ver2front" != "$2" ]; then
        [ "$ver1front" -gt "$ver2front" ] && return 11
        [ "$ver1front" -lt "$ver2front" ] && return 9

        [ "$ver1front" == "$1" ] || [ -z "$ver1back" ] && ver1back=0
        [ "$ver2front" == "$2" ] || [ -z "$ver2back" ] && ver2back=0
        do_version_check "$ver1back" "$ver2back"
        return $?
    else
        [ "$1" -gt "$2" ] && return 11 || return 9
    fi
}

# check version because 2.04 had a bug: it added a newline to the end
# of the file on every invocation

# version number last word of output, see
# http://stackoverflow.com/a/16620897/987623
astyle_version=$($ASTYLE --version | grep -oE '[^ ]+$')
minimum_astyle_version="2.5"
do_version_check "${astyle_version}" "${minimum_astyle_version}"
if [ $? -eq 9 ]; then
    echo "Found astyle version ${astyle_version}. Need at least ${minimum_astyle_version}!"
    exit 1
fi

RETURN=0
git diff --cached --name-status --diff-filter=ACMR |
    {
	    # Command grouping to workaround subshell issues. When the while loop is
	    # finished, the subshell copy is discarded, and the original variable
	    # RETURN of the parent hasn't changed properly.
	    while read STATUS FILE; do
            # regex matching
	        if [[ "$FILE" =~ \.(C|cxx|h)$ ]]; then
		        $ASTYLE $OPTIONS < $FILE > $FILE.beautified
		        md5sum -b $FILE | { read stdin; echo $stdin.beautified; } | md5sum -c --status -
		        if [ $? -ne 0 ]; then
			        echo "[commit rejected] $FILE does not respect the agreed coding standards." >&2
			        RETURN=1
		        fi
		        rm $FILE.beautified
	        fi
	    done

	    if [ $RETURN -eq 1 ]; then
		    echo "">&2
		    echo "Make sure to run 'make format' and review the changes *before* committing."  >&2
	    fi
	    exit $RETURN
    }
