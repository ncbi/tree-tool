#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Select rows by cut from a tsv-file"
  echo "#1: tsv-file name"
  echo "#2: comma-separated list of column numbers"
  echo "#3: make rows unique (0/1)"
  exit 1
fi
F=$1
C="$2"
UNIQ=$3


TMP=`mktemp`


head -1 $F | cut -f $C 
tail -n +2 $F | cut -f $C > $TMP
if [ $UNIQ == 1 ]; then
  sort -u $TMP
else
  cat $TMP
fi


rm $TMP*

