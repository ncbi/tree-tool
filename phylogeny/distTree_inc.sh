#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Build a distance tree incrementally"
  echo "Update: #1/"
  echo "Output: leaf_errors.{dm,txt}, tree.<DATE>, disagreement_nodes.txt, disagreement_nodes, gain_nodes, qual, core"
  echo "#1: incremental distance tree directory"
  echo "#2: seed (>=1)"
  echo "Time: O(n log^5(n))"
  exit 1
fi


if [ 1 == 1 ]; then  
# Time: O(n log^5(n))
while [ 1 == 1 ]; do
  if [ -e $1/stop ]; then
    echo "Stopped"
    exit 2
  fi
  
  ADD=`ls $1/new/ | wc -l`
  echo "# Add: $ADD  `date`  `date +%s`" >> $1/runlog  
  echo ""
  echo ""
	set +o errexit
  distTree_inc_new.sh $1 $2 
  S=$?
	set -o errexit
  if [ $S == 2 ]; then
    break
  fi
  if [ $S -ne 0 ]; then
    exit 1
  fi
done
  

echo ""
echo ""
echo "Complete optimization ..."
echo "# Complete optimization  `date`  `date +%s`" >> $1/runlog  
VER=`cat $1/version`
# Time: O(n log(n)) 
cp $1/tree $1/hist/tree.$VER
#
VER=$(( $VER + 1 ))
echo $VER > $1/version
# Time: O(n log^5(n))
makeDistTree  -threads 15  -data $1/  -variance lin \
  -optimize  -skip_len  -reinsert  \
  -output_tree $1/tree.new  -leaf_errors leaf_errors > $1/hist/makeDistTree.$VER
mv $1/tree.new $1/tree
tail -n +5 leaf_errors.dm | sort -k 2 -g -r > leaf_errors.txt
else
  VER=`cat $1/version`
fi 


distTree_inc_tree1_quality.sh $1


if [ -e $1/phen ]; then
	echo ""
	echo "Quality ..."
	tree_quality_phen.sh $1/tree $1/phen > $1/hist/tree_quality_phen.$VER 
	new_root=`grep '^New root: ' $1/hist/tree_quality_phen.$VER | sed 's/^New root: //1'`
	DATE=`date +%Y%m%d`

	echo ""
	echo ""
	echo "New root: $new_root"
	echo ""
	echo ""
	makeDistTree  -data $1/  -reroot_at "$new_root"  -output_tree tree.$DATE > /dev/null
	echo ""
	tree_quality_phen.sh tree.$DATE $1/phen > $1/hist/tree_quality_phen.$VER-rooted
	cat $1/hist/tree_quality_phen.$VER-rooted
	new_root=`grep '^New root: ' $1/hist/tree_quality_phen.$VER-rooted | sed 's/^New root: //1'`
	if [ "$new_root" ]; then
	  echo "Re-rooting must be idempotent"
	  exit 1
	fi
fi
