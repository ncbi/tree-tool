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


#set -x


# PAR
#CPP_DIR/phylogeny/database/combine_dissims.sh $REQ $GENOME "$FILE_NEW" $OUT 50 0.5 $INC/dissim_scale $INC/hmm-univ.stat 0 0.56 1  $LOG
 #                                                     1    2       3           4    5  6   7                 8                  9 10   11 12 


if [ $FILE_NEW ]; then
  error "New object processing is not implemented"
fi


TMP=`mktemp`


INC=`dirname $0`
GENOME_DIR=$INC/../genome


function req2file
{
  COL=$1  # 1|2
  #
  awk '{printf "'$GENOME_DIR'/%s/%s.prot-univ\n", $'$COL', $'$COL'};' $REQ
}


req2file 1 > $TMP.1
req2file 2 > $TMP.2
paste $TMP.1 $TMP.2 > $TMP.req

CPP_DIR/dissim/prot_collection2dissim  -log $LOG  -raw_power 0.57  $INC/hmm-univ.stat  $TMP.req $TMP.out

cut -f 1 $TMP.out | sed 's/\.prot-univ//1' | tr '/' '\t' | cut -f 5 > $TMP.out1
cut -f 2 $TMP.out | sed 's/\.prot-univ//1' | tr '/' '\t' | cut -f 5 > $TMP.out2
cut -f 3 $TMP.out | sed 's/^-nan$/inf/1' | sed 's/^nan$/inf/1'      > $TMP.out3

paste $TMP.out1 $TMP.out2 $TMP.out3 > $OUT 


rm -f $LOG
rm $TMP*

