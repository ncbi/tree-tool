#!/bin/bash
source CPP_DIR/bash_common.sh
if [ $# -ne 3 ]; then
  echo "Compute dissimilarities"
  echo "#1: file with pairs <Object1> <Object2>"
  echo "#2: output file"
  echo "#3: log (temporary)"
  exit 1
fi
REQ=$1
OUT=$2
LOG=$3


INC=`dirname $0`
# PAR
CPP_DIR/dissim/combine_dissims3.sh $REQ $OUT $INC/../genome 50 0.5 $INC/dissim_scale $INC/hmm-univ.stat 0 0.57 $LOG


