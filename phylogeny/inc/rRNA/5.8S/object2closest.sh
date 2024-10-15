#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find 100 closest sequences in the tree"
  echo "#1: Input sequence"
  echo "#2: sequence directory or ''"
  echo "#3: subset of objects (absolute pathname)"
  echo "#4: output file (absolute pathname)"
  exit 1
fi
NEW_OBJ=$1
DIR="$2"
OUT=$4


INC=`dirname $0`
if [ -z $DIR ]; then
  H=`file2hash $NEW_OBJ`
  DIR=$INC/../seq/$H
fi
CPP_DIR/genetics/kmerIndex_find $INC/seq.kmi $DIR/$NEW_OBJ 100 -qc > $OUT
