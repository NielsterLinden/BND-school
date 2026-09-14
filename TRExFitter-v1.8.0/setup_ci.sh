#!/bin/sh

function version_gt() { test "$(printf '%s\n' "$@" | sort -V | head -n 1)" != "$1"; }


if [ "${BASH_SOURCE[0]}" != "" ]; then
    # This should work in bash.
    _src=${BASH_SOURCE[0]}
elif [ "${ZSH_NAME}" != "" ]; then
    # And this in zsh.
    _src=${(%):-%x}
elif [ "${1}" != "" ]; then
    # If none of the above works, we take it from the command line.
    _src="${1/setup.sh/}/setup.sh"
else
    echo -e "\033[1;31mERROR:\033[0m Could not determine the base directory of TRExFitter, i.e. where \"setup.sh\" is located."
    echo -e "\033[1;31mERROR:\033[0m Can you give it to the source script as additional argument?"
    echo -e "\033[1;31mERROR:\033[0m For example: source ../setup.sh .."
    return 1
fi

location="$(cd -P "$(dirname "${_src}")" && pwd)"
unset _src

if [[ "$location" != "" ]]
then
  export PATH=$location:$PATH  # prepend to take priority over StatAnalysis-provided version
  # to be able to point to the config schema
  export TREXFITTER_HOME=$location
else
  export PATH=`pwd`:$PATH
  # to be able to point to the config schema
  export TREXFITTER_HOME=`pwd`
fi

export PATH=${TREXFITTER_HOME}/build/bin:${PATH}  # prepend to take priority over StatAnalysis
