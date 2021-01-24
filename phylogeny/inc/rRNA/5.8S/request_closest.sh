#!/bin/bash
source CPP_DIR/phylogeny/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Find 100 closest sequences in the tree"
  echo "#1: Input sequence"
  exit 1
fi
NEW_OBJ=$1


H=`file2hash $NEW_OBJ`

INC=`dirname $0`
CPP_DIR/genetics/kmerIndex_find $INC/seq.kmi $INC/../seq/$H/$NEW_OBJ 100 -qc
