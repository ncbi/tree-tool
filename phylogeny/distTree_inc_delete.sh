#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Delete a list of objects from an incremental distance tree"
  echo "Update: #1/"
  echo "Delete: #2"
  echo "#1: incremental distance tree directory"
  echo "#2: sorted and distinct list of objects to delete; to be moved into #1/hist/"
  exit 1
fi
INC=$1
DEL=$2


wc -l $DEL


# QC
sort -c -u $DEL


N=$( file2var $INC/threads 15 )
THREADS="-threads $N"

VER=$( cat $INC/version )

cp $INC/tree $INC/hist/tree.$VER
if [ $VER -gt 1 ]; then
  gzip $INC/hist/tree.$VER
fi

VER=$(( VER + 1 ))
echo $VER > $INC/version


section "Deleting from the tree"
VARIANCE=$( cat $INC/variance )
# Cf. distTree_inc_new.sh
$THIS/makeDistTree  $THREADS  -data $INC/  -variance $VARIANCE \
  -delete $DEL  \
  -optimize  -skip_len  -subgraph_iter_max 1 \
  -noqual \
  -output_dissim $INC/dissim.new \
  -output_tree $INC/tree.new > $INC/hist/makeDistTree-delete.$VER
  
mv $INC/tree.new $INC/tree

wc -l $INC/dissim
wc -l $INC/dissim.new
mv $INC/dissim.new $INC/dissim


section "indiscern"
cp /dev/null $INC/indiscern
$THIS/distTree_inc_dissim2indiscern.sh $INC $INC/dissim


section "Database"
$INC/objects_in_tree.sh $DEL null


section "$INC/new/"
$THIS/distTree_inc_new_cmd.sh $INC "rm -f" $DEL


mv $DEL $INC/hist/delete.$VER


section "QC"
$INC/qc.sh 0 


$THIS/distTree_inc_tree1_quality.sh $INC


