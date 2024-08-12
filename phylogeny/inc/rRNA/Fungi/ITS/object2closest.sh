#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find 100 closest sequences in the tree"
  echo "#1: new sequence"
  echo "#2: directory of #1 or ''"
  echo "#3: subset of sequence id's" 
  exit 1
fi
NEW_OBJ=$1
DIR="$2"
SUBSET=$3


INC=`dirname $0`
if [ -z $DIR ]; then
  DIR=$INC/../seq-long
fi
CPP_DIR/genetics/dna_closest.sh $DIR/$NEW_OBJ $INC/seq.fa "$SUBSET"




