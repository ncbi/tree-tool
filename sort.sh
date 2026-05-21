#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: File to remove trailing spaces and sort uniquely"
  exit 1
fi
F=$1


TMP=$( mktemp )
sed 's/^ *//1' $F | sed 's/ *$//1' | grep -vx "" | sort -u > $TMP
mv $TMP $1
