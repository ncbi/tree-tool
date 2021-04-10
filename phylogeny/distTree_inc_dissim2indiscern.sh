#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add indiscernibility pairs to #1/indiscern"
  echo "#1: incremental distance tree directory"
  echo "#2: dissimilarity triples: <obj1> <obj2> <dissim>"
  exit 1
fi
INC=$1
DISSIM=$2


awk '{if ($3 == 0)  print $1, $2};' $DISSIM >> $INC/indiscern
