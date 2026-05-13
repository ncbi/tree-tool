#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: File to sort (uniquely)"
  exit 1
fi

TMP=$( mktemp )
sort $1 -u > $TMP
mv $TMP $1
