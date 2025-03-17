#!/bin/bash
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Compute dissimilarities"
  echo "#1: file with pairs: <Object1> <Object2>"
  echo "#2: new object file or directory, or ''. Object name is basename"
  echo "#3: output file with triples:  <Object1> <Object2> <dissimilarity>"
  echo "#4: log (temporary)"
  exit 1
fi
REQ=$1
FILE_NEW="$2"
OUT=$3
LOG=$4


if [ -n "$FILE_NEW" ]; then
  error "New object $FILE_NEW"
fi


#TMP=$( mktemp )
#echo $TMP
#set -x


INC=$THIS
CPP_DIR/trav $REQ "CPP_DIR/dissim/seq2dissim $INC/../seq/%1 $INC/../seq/%2  -prot_name aa  -global -blosum62 -log $LOG | grep -w '^min_edit:' | cut -f 2 | sed 's/^/%f\t/1'" > $OUT
rm -f $LOG


#rm $TMP*

