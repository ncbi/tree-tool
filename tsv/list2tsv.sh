#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Convert a list into a tsv-file"
  echo "#1: list file"
  echo "#2: column name"
  exit 1
fi
F=$1
COL="$2"


echo "#"$COL 
sort -u $F 
