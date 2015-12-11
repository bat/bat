# get directory of this script when sourcing from somewhere else
# http://stackoverflow.com/a/179231/987623
SCRIPT_PATH="${BASH_SOURCE[0]}";
if ([ -h "${SCRIPT_PATH}" ]) then
  while([ -h "${SCRIPT_PATH}" ]) do SCRIPT_PATH=`readlink "${SCRIPT_PATH}"`; done
fi
pushd . > /dev/null
cd `dirname ${SCRIPT_PATH}` > /dev/null
SCRIPT_PATH=`pwd`;
popd  > /dev/null

export LD_LIBRARY_PATH=${SCRIPT_PATH}:$LD_LIBRARY_PATH
