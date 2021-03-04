#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then 
  echo "Create a .dm-file with two-way attributes"
  echo "#1: prots_pair2dissim output"
  echo "#2: hmm-univ.list"
  echo "#3: Output .dm-file"
  exit 1
fi
IN=$1
HMMLIST=$2
OUT=$3


TMP=`mktemp`


OBJNUM=`cat $IN | wc -l`
echo "OBJNUM $OBJNUM NAME NOMULT" > $OUT
echo "ATTRIBUTES" >> $OUT
cat $HMMLIST | sed 's/$/ Positive2 6' >> $OUT
echo "DATA" >> $OUT

cut -f 1 -d ' ' $IN | tr '-' '\t' > $TMP.pairs
cut -f 1 $TMP.pairs > $TMP.list
cut -f 2 $TMP.pairs > $TMP.list
$THIS/../uniq.sh $TMP.list
cat $TMP.list >> $OUT

PAIRS=`cat $IN | wc -l`
ATTRS=`cat $HMMLIST | wc -l`
echo "PAIR_DATA $PAIRS ATTRS" >> $OUT
cat $HMMLIST | tr '\n' ' ' >> $OUT
echo "" >> $OUT

sed 's/-/ /1' $IN >> $OUT


rm $TMP*
