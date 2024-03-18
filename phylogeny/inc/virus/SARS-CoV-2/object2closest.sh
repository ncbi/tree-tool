#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find closest genomes among #3"
  echo "#1: Virus.id"
  echo "#2: new object directory or ''"
  echo "#3: subset of Genome.id's"
  exit 1
fi
NEW_OBJ=$1
DIR="$2"
SUBSET=$3


#set -o xtrace 


INC=`dirname $0`

if [ -z $DIR ]; then
  H=`CPP_DIR/file2hash $NEW_OBJ`
  DIR=$INC/../mut.dna/$H
fi


TMP=`mktemp`


SIZE=30000  # = ln(tree size)  * 2000
CPP_DIR/index_find $INC/../mut.dna $INC/../mut.index $SIZE $DIR/$NEW_OBJ -large > $TMP
if [ -s $TMP ]; then
  grep -f $INC/../deleted.all -vx $TMP > $TMP.out
  head -100 $TMP.out
else
  if true; then  
    REAL_REF=`readlink $INC/../ref`
    basename $REAL_REF
  else
    echo "MW737421.1"
  fi
fi


rm $TMP*

