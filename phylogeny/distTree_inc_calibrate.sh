#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "For #1.list compute dissimilarities #1.dm, distance tree #1.tree and evaluate it"
  echo "#1: list of objects (file prefix)"
  echo "#2: incremental distance tree directory (no updating)"
  exit 1
fi
F=$1
INC=$2


LIST=$F.list
if [ ! -e $LIST ]; then
  error "File $LIST does not exist"
fi

if [ ! -e $INC/phen ]; then
  error "Directory $INC/phen does not exist"
fi



section "Dissimilarities"
TMP=`mktemp`
echo $TMP

$THIS/../list2pairs $LIST > $TMP.req
$THIS/distTree_inc_request2dissim.sh $INC $TMP.req $TMP.dissim
$THIS/../dm/pairs2dm $TMP.dissim 1 "cons" 6  -distance > $F.dm

rm $TMP*


section "Tree"
N=15
if [ -e $INC/threads ]; then
  N=`cat $INC/threads`
fi
THREADS="-threads $N"

VARIANCE=`cat $INC/variance`

$THIS/makeDistTree  $THREADS  -data $F  -dissim_attr "cons"  -variance $VARIANCE  -optimize  -subgraph_iter_max 5  -noqual  -output_tree $F.tree


section "Quality"
LARGE=0
if [ -e $INC/large ]; then
  LARGE=1
fi
$THIS/tree_quality_phen.sh $F.tree "" $INC/phen $LARGE 1 ""

