#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Sort a tsv-file"
  echo "#1: tsv-file name"
  echo "#2: sort parameters"
  exit 1
fi
F=$1
P="$2"


if head -1 $F | grep -s '^#'; then
  grep    '^#' $F
  grep -v '^#' $F | sort -t $'\t' $P
else
  printf '#'
  head -1 $F
  tail -n +2 $F | sort -t $'\t' $P
fi


