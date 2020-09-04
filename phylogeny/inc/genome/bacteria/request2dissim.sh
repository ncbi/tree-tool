#!/bin/bash
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Compute dissimilarities for pairs of objects"
  echo "#1: input dissimilarity requests (pairs of objects)"
  echo "#2: new object file or directory, or ''. Object name is basename"
  echo "#3: output dissimilarities added to the pairs of objects"
  echo "#4: error log"
  exit 1
fi
REQ=$1
FILE_NEW="$2"
OUT=$3
LOG=$4


if [ -n "$FILE_NEW" ]; then
  error "New object is not implemented"
fi

INC=`dirname $0`
# PAR
CPP_DIR/dissim/combine_dissims3.sh $REQ $OUT $INC/../genome 50 0.5 $INC/dissim_scale $INC/hmm-univ.stat 0 0.57 $LOG


