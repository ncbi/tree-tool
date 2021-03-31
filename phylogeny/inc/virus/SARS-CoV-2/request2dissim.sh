#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "$0"
  echo "#1: dissimilarity requests (input)"
  echo "#2: new object or ''"
  echo "#3: dissim (output)"
  echo "#4: log"
  exit 1
fi
REQUEST=$1
NEW_OBJ=$2
DISSIM=$3
LOG=$4


if [ ! -s $REQUEST ]; then
  error "Empty request" >> $LOG
fi

if [ -n "$NEW_OBJ" ]; then
  error "NEW_OBJ is not implemented"
fi

INC=`dirname $0`
#dna_pairs2dissim  -log $LOG  -diff  -band 100  $REQUEST $INC/../seq 20000 $DISSIM
CPP_DIR/dissim/feature_request2dissim $REQUEST $INC/../mut.feature  -optional_weight 0.2  -virus  -log $LOG  > $DISSIM

if [ ! -s $DISSIM ]; then
  error "Empty dissim" >> $LOG
fi

rm -f $LOG
