#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: incremental distance tree directory"
  exit 1
fi
INC=$1


TMP=`mktemp`
#echo $TMP


awk '{if ($3 == 0 && $1 < $2)  print $1, $2};' $INC/dissim >  $TMP
awk '{if ($3 == 0 && $1 > $2)  print $2, $1};' $INC/dissim >> $TMP
sort -u $TMP > $TMP.dissim
#wc -l $TMP.dissim

awk '{if ($1 < $2)  print $1, $2};' $INC/indiscern >  $TMP
awk '{if ($1 > $2)  print $2, $1};' $INC/indiscern >> $TMP
sort -u $TMP > $TMP.indiscern
#wc -l $TMP.indiscern

diff $TMP.dissim $TMP.indiscern


rm $TMP*
