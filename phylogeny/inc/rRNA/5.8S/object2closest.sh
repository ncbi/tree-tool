#!/bin/bash --noprofile
source $PANFS/code/cpp/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find 100 closest sequences in the tree"
  echo "#1: Input sequence"
  echo "#2: sequence directory or ''"
  exit 1
fi
NEW_OBJ=$1
DIR="$2"


INC=`dirname $0`
if [ -z $DIR ]; then
  H=`file2hash $NEW_OBJ`
  DIR=$INC/../seq/$H
fi
$PANFS/code/cpp/genetics/kmerIndex_find $INC/seq.kmi $DIR/$NEW_OBJ 100 -qc
