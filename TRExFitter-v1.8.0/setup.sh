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


# Figure out which version of Linux OS is being used from OS CPE name
cpe_name=`cat /etc/os-release | grep CPE_NAME`
cpe_name=${cpe_name:9}

if [[ $cpe_name == *"linux:9"* ]]; then
    os_version=9
elif  [[ $cpe_name == *"centos:7"* ]]; then
    os_version=7
else
    echo "Cannot recognise your OS, will try EL9 setup..."
    os_version=9
fi

# Setup ROOT and gcc
# added back by Michele
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet

asetup StatAnalysis,0.7.2  # ROOT v6-38-00

if [ "${ROOTSYS}" = "" ]; then
   echo -e "\033[41;1;37m Error initializing ROOT. ROOT is not set up. Please check. \033[0m"
else
   echo -e "\033[42;1;37m ROOT has been set to: *${ROOTSYS}* \033[0m"
fi

alias macro="root -l -b -q"

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

if [ ! -f $TREXFITTER_HOME/logo.txt ]; then
  echo -e "\033[1;31mWARNING:\033[0m \$TREXFITTER_HOME environmental variable not set properly"
  echo -e "\033[1;31mWARNING:\033[0m call this script with the path to the TRExFitter directory as an additional argument"
fi

# Check if the xrootfit code exists
if [ ! "$(ls -A ${TREXFITTER_HOME}/xroofit)" ]; then
  echo -e "\033[1;31mERROR:\033[0m xroofit directory does not exist or is empty. "
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule init' (first time use) in the base directory of TRExFitter"
  echo -e "\033[1;31mERROR:\033[0m Followed with 'git submodule update' in the base directory of TRExFitter"
  return
fi

# Check if the Blinder code exists
if [ ! "$(ls -A ${TREXFITTER_HOME}/Blinder)" ]; then
  echo -e "\033[1;31mERROR:\033[0m Blinder directory does not exist or is empty. "
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule init' (first time use) in the base directory of TRExFitter"
  echo -e "\033[1;31mERROR:\033[0m Followed with 'git submodule update' in the base directory of TRExFitter"
  return
fi

# Check if the CommomSmoothing code exists
if [ ! "$(ls -A ${TREXFITTER_HOME}/SystematicSmoothingTool)" ]; then
  echo -e "\033[1;31mERROR:\033[0m SystematicSmoothingTool directory does not exist or is empty. "
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule init' (first time use) in the base directory of TRExFitter"
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule update' in the base directory of TRExFitter"
  return
fi

# Check if the UnfoldingCode code exists
if [ ! "$(ls -A ${TREXFITTER_HOME}/UnfoldingCode)" ]; then
  echo -e "\033[1;31mERROR:\033[0m UnfoldingCode directory does not exist or is empty. "
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule init' (first time use) in the base directory of TRExFitter"
  echo -e "\033[1;31mERROR:\033[0m Followed with 'git submodule update' in the base directory of TRExFitter"
  return
fi

# Check if the yaml-cpp code exists
if [ ! "$(ls -A ${TREXFITTER_HOME}/yaml-cpp)" ]; then
  echo -e "\033[1;31mERROR:\033[0m yaml-cpp directory does not exist or is empty. "
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule init' (first time use) in the base directory of TRExFitter"
  echo -e "\033[1;31mERROR:\033[0m Followed with 'git submodule update' in the base directory of TRExFitter"
  return
fi

echo "Setting up cmake with: lsetup cmake"
lsetup cmake

export PATH=${TREXFITTER_HOME}/build/bin:${PATH}  # prepend to take priority over StatAnalysis
alias trex-make='cd build/; make -j4; cd ../'
alias trex-clean='rm -rf build/'
alias trex-build='mkdir build && cd build; cmake ../; make -j4; cd ../'

version=`cat ${TREXFITTER_HOME}/version.txt`
echo -e "\n\e[1m${version} ready (if everything went smoothly)\e[0m"
echo -e "************************************************************************"
echo -e "First time?"
echo -e "  \e[1mmkdir -p build && cd build && cmake ../ && cmake --build ./ && cd ..\e[0m"
echo -e "To recompile:"
echo -e "  \e[1mtrex-make\e[0m"
echo -e "To run the code:"
echo -e "  \e[1mtrex-fitter <action(s)> <config file> [<options>]\e[0m"
echo -e "See README.md for more details."
echo -e "************************************************************************"
echo -e "\n"
