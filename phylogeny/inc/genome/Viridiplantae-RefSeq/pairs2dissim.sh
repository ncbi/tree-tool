#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Compute dissimilarities for pairs of objects"
  echo "#1: input dissimilarity requests (pairs of objects)"
  echo "#2: new object path or ''. Object name is basename"
  echo "#3: output dissimilarities added to the pairs of objects"
  echo "#4: error log"
  exit 1
fi
REQ=$1
FILE_NEW="$2"
OUT=$3
LOG=$4


#set -x


if [ $FILE_NEW ]; then
  exit 2
fi


INC=`dirname $0`
GENOME=$INC/../genome
# PAR
#CPP_DIR/phylogeny/database/combine_dissims.sh $REQ $GENOME "$FILE_NEW" $OUT 200 0.1 $INC/dissim_scale $INC/hmm-univ.stat 1 0.38 1  $LOG  
#                                                  1    2       3           4    5   6   7                 8                  9 10   11 12 

TMP=`mktemp`

awk '{printf "'$GENOME'/%s/%s.prot-univ '$GENOME'/%s/%s.prot-univ\n", $1, $1, $2, $2};' $REQ > $TMP
CPP_DIR/dissim/prot_collection2dissim  -log $LOG  -raw_power 0.38  -blosum62  -coeff 1.0  $INC/hmm-univ.stat $TMP $TMP.res
cut -f 3 $TMP.res > $TMP.dissim
paste $REQ $TMP.dissim > $OUT

rm $TMP*


rm -f $LOG

