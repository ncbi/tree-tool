#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# != 1 ]; then
  echo "#1: File to sort and uniq"
  exit 1
fi

TMP=$( mktemp )
< $1 tr '\t' ' ' | sed 's/ *$//1' | sed 's/  / /g' | sort -u | tr ' ' '\t' > $TMP
mv $TMP $1

