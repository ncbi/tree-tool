#! /bin/echo Please-source

set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C


export NOCOLOR='\033[0m'
export RED='\033[1;31m'
export GREEN='\033[1;32m'
export BROWN='\033[33m'
export YELLOW='\033[1;33m'
export BLUE='\033[1;34m'
export MAGENTA='\033[35m'


function error
{
	local MSG="$1"
  echo ""
	echo -e "${RED}ERROR${NOCOLOR}\07"
	echo -e "$MSG"
	exit 1
}


function section
{
	local MSG="$1"
  echo ""
  echo -e "${GREEN}${MSG}${NOCOLOR}"
}


function super_section
{
  local MSG="$1"
  local S="*** $MSG ***"
  local T=${S//?/*}  
  echo ""
  section "$T\n$S\n$T"
  date
}


function warning
{
	local MSG="$1"
  echo -e "${YELLOW}${MSG}${NOCOLOR}"
}


function success
{
  echo ""
  echo ""
  echo -e "${BLUE}SUCCESS${NOCOLOR}"
  if [ $SHLVL -le 2 ]; then
    printf "\07\07"
    sleep 1
    printf "\07\07"
    sleep 1
    printf "\07\07"
  fi
  date
 #exit 0
}


function comment
{
	local MSG="$1"
  echo -e "${BROWN}${MSG}${NOCOLOR}" >> /dev/stderr  # MAGENTA
}


function file2var
{
  local F=$1
  local DEFAULT="$2"
  if [ -s $F ]; then
    local R=$( cat $F )
    if [ "$R" ]; then
      echo "$R"
      return
    fi
  fi
  echo $DEFAULT
}


function differ 
{
  local F1=$1
  local F2=$2
  ( set -x && diff $F1 $F2 )
}


