#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Select rows by awk from a .tsv-file"
  echo "#1: .tsv-file name"
  echo "#2: file with column names without quotes"
  echo "#3: make rows unique (0/1) and non-empty"
  exit 1
fi
F=$1
NAMES=$2
UNIQ=$3


C=$( $THIS/../trav $NAMES -noprogress "$THIS/tsv_colname2num.sh $F %Q%f%Q" | tr '\n' ',' | sed 's/,$//1' )
$THIS/tsv_cut.sh $F $C $UNIQ

