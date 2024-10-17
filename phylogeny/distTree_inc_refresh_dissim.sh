#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Update #1/dissim"
  echo "#1: incremental distance tree directory"
  echo "#2: compute dissimilarity requests (0/1)"
  echo "#3: save old #1/dissim in #1/dissim.<version>.gz (0/1)"
  exit 1
fi
INC=$1
REQ=$2
SAVE=$3


VER=$( cat $INC/version )
THREADS=$( file2var $INC/threads 15 )


wc -l $INC/dissim


if [ $SAVE == 1 ]; then
  section "Saving $INC/dissim to $INC/dissim.$VER"
  if [ -e $INC/dissim.$VER ]; then
    error "$INC/dissim.$VER exists"
  fi
  cp $INC/dissim $INC/dissim.$VER && gzip $INC/dissim.$VER &
fi


TMP=$( mktemp )
comment $TMP


wc -l $INC/dissim

if [ $REQ == 1 ]; then
  section "Computing dissimilarity requests"
  $THIS/distTree_refresh_dissim $INC $TMP.req $INC/dissim  -threads $THREADS
else
  cut -f 1,2 $INC/dissim > $TMP.req
fi

section "Computing dissimilarities"
$THIS/distTree_inc_request2dissim.sh $INC $TMP.req $TMP.dissim-add

section "Updating $INC/dissim and $INC/indiscern"
if [ $REQ == 1 ]; then
  $THIS/distTree_inc_dissim2indiscern.sh $INC $TMP.dissim-add
  cat $TMP.dissim-add >> $INC/dissim
else
  mv $TMP.dissim-add $INC/dissim
  cp /dev/null $INC/indiscern
  $THIS/distTree_inc_dissim2indiscern.sh $INC $INC/dissim
  section "Updating $INC/tree"
  cp $INC/tree $INC/hist/tree.$VER
  gzip $INC/hist/tree.$VER
  VARIANCE=$( cat $INC/variance )
  $THIS/makeDistTree  -data $INC/  -output_tree $TMP.tree  -threads $THREADS  -variance $VARIANCE  -fix_discernible  -optimize  -skip_len  -subgraph_iter_max 2 > $INC/hist/makeDistTree.$VER
  mv $TMP.tree $INC/tree
fi

section "QC $INC/indiscern"
$THIS/distTree_inc_indiscern_qc.sh $INC 0 0

VER=$(( VER + 1 ))
echo $VER > $INC/version


rm $TMP*
wait
