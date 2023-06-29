#!/bin/bash --noprofile
source $PANFS/code/cpp/bash_common.sh
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


NEW=""
if [ -n "$FILE_NEW" ]; then
  NAME=`head -1 $FILE_NEW | sed 's/^>//1' | cut -f 1 -d ' '`
  NEW="-name_new $NAME  -file_new $FILE_NEW"
fi

INC=`dirname $0`
$PANFS/code/cpp/dissim/dna_pair2dissim  -log $LOG  -coeff 0.0066887  $NEW  $REQUEST $INC/../seq 300 $DISSIM
  # was: 0.00634 
  
rm -f $LOG


