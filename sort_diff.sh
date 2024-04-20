#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "sort and diff two text files"
  echo "#1: text file 1 (unique lines)"
  echo "#2: text file 2 (unique lines)"
  exit 1
fi
F1=$1
F2=$2


TMP=$( mktemp )

sort $F1 > $TMP.1
sort -cu $TMP.1

sort $F2 > $TMP.2
sort -cu $TMP.2

section "$F1 minus $F2:"
$THIS/setMinus $TMP.1 $TMP.2

section "$F2 minus $F1:"
$THIS/setMinus $TMP.2 $TMP.1


rm $TMP*
