#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Sort s tsv-file"
  echo "#1: tsv-file name"
  echo "#2: sort parameters"
  exit 1
fi
F=$1
P=$2


head -1 $F
tail -n +2 $F | sort -t $'\t' $P
