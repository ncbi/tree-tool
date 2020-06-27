#!/bin/bash
source CPP_DIR/bash_common.sh
if [ $# -ne 3 ]; then
  echo "Compute dissimilarities"
  echo "#1: file with pairs <Object1> <Object2>"
  echo "#2: Output file"
  echo "#3: log (temporary)"
  exit 1
fi
REQ=$1
OUT=$2
LOG=$3


INC=`dirname $0`
CPP_DIR/dissim/combine_dissims.sh $REQ $OUT $INC/../genome  10 0.1 $INC/dissim_scale $INC/hmm-univ.stat 1 0.57 $LOG
#                                                 1    2    3               4  5   6                 7                  8 9    10




