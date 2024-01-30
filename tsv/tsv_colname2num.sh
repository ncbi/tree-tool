#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "tsv_cat for a directory of .tsv-files"
  echo "#1: .tsv-file"
  echo "#2: column name"
  exit 1
fi
F=$1
COL="$2"


$THIS/tsv_schema $F | cut -f 2 | grep -xn "$COL" | tr ':' '\t' | cut -f 1

