#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Select rows by cut from a tsv-file"
  echo "#1: tsv-file name"
  echo "#2: comma-separated list of column numbers"
  exit 1
fi
F=$1
C="$2"


head -1 $F | cut -f $C 
tail -n +2 $F | cut -f $C | sort -u
