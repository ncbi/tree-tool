#!/bin/bash --noprofile
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


INC=`dirname $0`
# PAR
CPP_DIR/dissim/combine_dissims.sh $REQ $INC/../genome "$FILE_NEW" $OUT 10 0.1 $INC/dissim_scale $INC/hmm-univ.stat 1 0.65 1  1  $LOG
#                                         1    2              3           4    5  6   7                  8                 9 10   11 12 13
