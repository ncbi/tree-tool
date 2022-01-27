#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Select rows by awk from a tsv-file"
  echo "#1: tsv-file name"
  echo "#2: awk selection parameters"
  exit 1
fi
F=$1
Q=\'
P="$2"


head -1 $F
tail -n +2 $F | awk -F '\t' "$P"
