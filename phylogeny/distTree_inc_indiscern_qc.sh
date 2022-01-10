#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "#1: incremental distance tree directory"
  echo "#2: compute dissimilarities (0/1)"
  echo "#3: verbose (0/1)"
  exit 1
fi
INC=$1
COMP=$2
VERB=$3


TMP=`mktemp`
if [ $VERB == 1 ]; then
  echo $TMP
  set -x
fi


awk '{if ($3 == 0 && $1 < $2)  print $1, $2};' $INC/dissim >  $TMP
awk '{if ($3 == 0 && $1 > $2)  print $2, $1};' $INC/dissim >> $TMP
sort -u $TMP > $TMP.dissim
#wc -l $TMP.dissim

awk '{if ($1 < $2)  print $1, $2};' $INC/indiscern >  $TMP
awk '{if ($1 > $2)  print $2, $1};' $INC/indiscern >> $TMP
sort -u $TMP > $TMP.indiscern
#wc -l $TMP.indiscern

diff $TMP.dissim $TMP.indiscern


$THIS/tree2indiscern $INC/tree > $TMP.tree2indiscern
$THIS/../connectPairs $TMP.tree2indiscern $TMP.1  -pairs 
$THIS/../sort.sh $TMP.1

$THIS/tree2obj.sh $INC/tree > $TMP.list
$THIS/../connectPairs $INC/indiscern $TMP.2  -pairs  -subset $TMP.list
$THIS/../sort.sh $TMP.2

diff $TMP.1 $TMP.2


if [ $COMP == 1 ]; then
  $THIS/distTree_inc_request2dissim.sh $INC $INC/indiscern $TMP.dissim
  awk '$3 != 0' $TMP.dissim > $TMP.bad
  if [ -s $TMP.bad ]; then
    wc -l $TMP.bad
    exit 1
  fi
fi


rm $TMP*
