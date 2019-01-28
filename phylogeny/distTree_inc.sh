#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Build a distance tree incrementally"
  echo "Update: #1/"
  echo "Output: leaf_errors.{dm,txt}, tree.<DATE>, disagreement_nodes.txt, disagreement_nodes, gain_nodes, qual, core"
  echo "#1: incremental distance tree directory"
  echo "#2: new.list | ''"
  echo "Time: O(n log^5(n))"
  exit 1
fi
INC=$1
NEW=$2


if [ -e $NEW ]; then
  wc -l $NEW
  $THIS/../trav $NEW "cp /dev/null $INC/new/%f"
fi


VARIANCE=`cat $INC/variance`


if [ 1 == 1 ]; then   
# Time: O(n log^5(n))
while [ 1 == 1 ]; do
  if [ -e $INC/stop ]; then
    echo ""
    echo '*** STOPPED ***'
    exit 2
  fi
  
  if [ -e $INC/skip ]; then
    echo ""
    echo '*** SKIPPED ***'
    break
  fi
  
  ADD=`ls $INC/new/ | wc -l`
  echo "# Add: $ADD  `date`  `date +%s`" >> $INC/runlog  
  echo ""
  echo ""
	set +o errexit
  $THIS/distTree_inc_new.sh $INC  
  S=$?
	set -o errexit
  if [ $S == 2 ]; then
    break
  fi
  if [ $S -ne 0 ]; then
    exit $S
  fi
done
  

echo ""
echo ""
echo "Complete optimization ..."
echo "# Complete optimization  `date`  `date +%s`" >> $INC/runlog  
VER=`cat $INC/version`
# Time: O(n log(n)) 
cp $INC/tree $INC/hist/tree.$VER
gzip $INC/hist/tree.$VER
#
VER=$(( $VER + 1 ))
echo $VER > $INC/version
# Time: O(n log^5(n))
# PAR
$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE \
  -optimize  -skip_len  -reinsert  -subgraph_iter_max 10 \
  -output_tree $INC/tree.new  -leaf_errors leaf_errors > $INC/hist/makeDistTree-complete-inc.$VER
mv $INC/tree.new $INC/tree
tail -n +5 leaf_errors.dm | sort -k 2 -g -r > leaf_errors.txt

echo ""
$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE  -qc  -noqual > $INC/hist/makeDistTree-qc.$VER
else
  VER=`cat $INC/version`
fi 


$THIS/distTree_inc_tree1_quality.sh $INC



if [ -e $INC/phen ]; then
	DATE=`date +%Y%m%d`

	echo ""
	echo "Root and quality ..."
	$THIS/tree_quality_phen.sh $INC/tree "" $INC/phen > $INC/hist/tree_quality_phen.$VER 
	cat $INC/hist/tree_quality_phen.$VER 
	OLD_ROOT=`grep '^Old root: ' $INC/hist/tree_quality_phen.$VER | sed 's/^Old root: //1'`
	NEW_ROOT=`grep '^New root: ' $INC/hist/tree_quality_phen.$VER | sed 's/^New root: //1'`

	echo ""
	echo "Setting root and sorting ..."
  if [ ! "$NEW_ROOT" ]; then
    NEW_ROOT=$OLD_ROOT
  fi
	$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE  -reroot_at "$NEW_ROOT"  -noqual  -output_tree tree.$DATE > /dev/null  	
	
	echo ""
	echo "Names ..."
	$THIS/tree2names.sh tree.$DATE $INC/phen > $INC/hist/tree2names.$VER
fi


echo ""
date

