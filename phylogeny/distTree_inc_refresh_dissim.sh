#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Compute new #1/dissim given #1/tree; save old #1/dissim in #1/dissim.<version>"
  echo "#1: incremental distance tree directory"
  exit 1
fi
INC=$1


VER=`cat $INC/version`

N=15
if [ -e $INC/threads ]; then
  N=`cat $INC/threads`
fi
THREADS="-threads $N"


section "Saving $INC/dissim to $INC/dissim.$VER"
if [ -e $INC/dissim.$VER ]; then
  error "$INC/dissim.$VER exists"
fi
wc -l $INC/dissim
cp $INC/dissim $INC/dissim.$VER
gzip $INC/dissim.$VER &


TMP=`mktemp`
echo $TMP


section "Computing dissimilarity requests"
$THIS/distTree_refresh_dissim $INC $TMP.req $INC/dissim $THREADS
wc -l $INC/dissim

section "Computing dissimilarities"
$THIS/distTree_inc_request2dissim.sh $INC $TMP.req $TMP.dissim-add

section "Updating $INC/indiscern"
$THIS/distTree_inc_dissim2indiscern.sh $INC $TMP.dissim-add
cat $TMP.dissim-add >> $INC/dissim


rm $TMP*
wait
