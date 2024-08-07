#!/bin/bash --noprofile
source $PANFS/code/cpp/bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find 100 closest sequences in the tree"
  echo "#1: new sequence"
  echo "#2: directory of #1 or ''"
  echo "#3: subset of sequence id's"  # ignored ??
  exit 1
fi
NEW_OBJ=$1
DIR="$2"


INC=`dirname $0`
if [ -z $DIR ]; then
  DIR=$INC/../seq
fi
$PANFS/code/cpp/genetics/dna_closest.sh $DIR/$NEW_OBJ $INC/seq.fa


