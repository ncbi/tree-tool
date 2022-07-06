#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add sequential numbers to the end of a tsv-file"
  echo "#1: tsv-file name"
  echo "32: sequence column name"
  exit 1
fi
F=$1
COL="$2"


head -1 $F | sed "s/$/\t$COL/1"
tail -n +2 $F | awk -F '\t' '{OFS="\t"; print $0, NR};'
