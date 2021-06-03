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
REQUEST=$1
FILE_NEW="$2"
DISSIM=$3
LOG=$4


if [ $FILE_NEW ]; then
  error "$0 is not implemented"
fi


TMP=`mktemp`


INC=`dirname $0`
awk '{printf "'$INC'/../genome/%s/%s.prot-univ '$INC'/../genome/%s/%s.prot-univ\n", $1, $1, $2, $2};' $REQUEST > $TMP
CPP_DIR/dissim/prot_collection2dissim  $INC/hmm-univ.stat $TMP $DISSIM  -raw_power 0.6  -blosum62  -coeff 2.25  -log $LOG


rm -f $LOG
rm $TMP*
