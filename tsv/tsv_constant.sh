#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Add a constant column to a tsv-file"
  echo "#1: tsv-file name"
  echo "#2: column name"
  echo "#3: value"
  exit 1
fi
F=$1
COL="$2"
VAL="$3"


head -1 $F | sed "s/$/\t$COL/1"
tail -n +2 $F | awk -F '\t' '{OFS="\t"; print $0, "'$VAL'"};'
