#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Print #2 closest objects"
  echo "#1: incremental distance tree"
  echo "#2: target object"
  echo "#3: number of closest objects to #2"
  echo "#4: print dissimilaritiies? (0/1)"
  exit 1
fi
INC=$1
OBJ=$2
NUM=$3
FULL=$4


TMP=$( mktemp )


grep -w $OBJ $INC/dissim | sort -k 3 -g > $TMP
head -$NUM $TMP > $TMP.head
if [ $FULL == 1 ]; then
  cat $TMP.head
else
  < $TMP.head tr ' ' '\t' | cut -f 1,2 | tr '\t' '\n' | grep -vw $OBJ
fi


rm $TMP*
