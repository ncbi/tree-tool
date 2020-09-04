#!/bin/bash
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Find 100 closest sequences in the tree"
  echo "#1: Input sequence"
  exit 1
fi
NEW_OBJ=$1


INC=`dirname $0`
CPP_DIR/genetics/dna_closest.sh $INC/../seq/$NEW_OBJ $INC/seq.fa




