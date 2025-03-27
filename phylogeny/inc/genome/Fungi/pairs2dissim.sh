#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Compute dissimilarities for pairs of objects"
  echo "#1: input dissimilarity requests (pairs of objects)"
  echo "#2: new object path or ''. Object name is basename"
  echo "#3: output dissimilarities added to the pairs of objects (absolute pathname)"
  echo "#4: error log"
  exit 1
fi
REQ=$1
FILE_NEW="$2"
OUT=$3
LOG=$4


#set -x


INC=$( dirname $0 )
GENOME=$INC/../genome
# PAR
CPP_DIR/phylogeny/database/combine_dissims.sh $REQ $GENOME "$FILE_NEW" $OUT 200 0.1 $INC/dissim_scale $INC/hmm-univ.stat 1 0.80 1  $LOG
#                                                 1    2       3           4    5   6   7                 8                  9 10   11 12 
# was: #5 = 10
#      #9 = 1

#CPP_DIR/trav  -step 1  $REQ "CPP_DIR/phylogeny/database/symbet.sh $GENOME %f 20 2> /dev/null | sed 's/^/%f\t/1'"  -log $LOG > $OUT


if false; then
  TMP=$( mktemp )
  #echo $TMP 
  function req2file
  {
    local REQ_=$1
    local COL=$2  # 1|2
    local SUF=$3
    #
    cut -f $COL $REQ_ > $TMP.req2file_req
    CPP_DIR/file2hash $TMP.req2file_req -file -append  -log $LOG  -noprogress | awk '{printf "'$GENOME'/%s/%s/%s.'$SUF'\n", $1, $2, $2};'
  }
  req2file $REQ 1 "prot-univ" > $TMP.f1
  req2file $REQ 2 "prot-univ" > $TMP.f2
  paste $TMP.f1 $TMP.f2 > $TMP.req-univ
  CPP_DIR/dissim/prot_collection2dissim  -log $LOG  -power 0.57  $INC/hmm-univ.stat  $TMP.req-univ $TMP.univ
  cut -f 3 $TMP.univ > $TMP.dissim-univ_1
  paste $REQ $TMP.dissim-univ_1 > $OUT
  rm $TMP*
fi


rm -f $LOG

