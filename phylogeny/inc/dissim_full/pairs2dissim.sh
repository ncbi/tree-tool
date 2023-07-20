#!/bin/bash
THIS=`dirname $0`
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


TMP=`mktemp`
#echo $TMP
#set -x


awk '{if ($1 < $2)  printf "%s\t%s\n", $1, $2};' $REQ >  $TMP
awk '{if ($1 > $2)  printf "%s\t%s\n", $2, $1};' $REQ >> $TMP
sort -u $TMP | sed 's/\t/-/1' > $TMP.req
join  -1 1  -2 1  $TMP.req $THIS/../dissim_full > $TMP.out
cut  -f 1  -d ' '  $TMP.out > $TMP.found
CPP_DIR/setMinus $TMP.req $TMP.found | sed 's/$/ inf/1' >> $TMP.out

sed 's/-/ /1' $TMP.out | tr ' ' '\t' > $OUT


rm $TMP*

