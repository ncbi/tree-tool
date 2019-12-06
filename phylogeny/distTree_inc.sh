#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Build a distance tree incrementally"
  echo "Update: #1/"
  echo "Output: leaf_errors.{dm,txt}, tree.<DATE>, disagreement_nodes[.txt], disagreement_objects, gain_nodes, qual, qual.raw"
  echo "#1: incremental distance tree directory"
  echo "#2: add new (0 - no, 1 - almost all, 2 - all)"
  echo "Time: O(n log^4(n))"
  exit 1
fi
INC=$1
NEW_PAR=$2


echo "QC ..."
$INC/qc.sh go
echo ""


if [ 1 == 1 ]; then  
if [ $NEW_PAR -gt 0 ]; then
  # Time: O(n log^4(n))
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
    
    NEW=`ls $INC/new/ | wc -l`
    VER=`cat $INC/version`
    echo "$VER  # New: $NEW  `date`  `date +%s`" >> $INC/runlog  
    echo ""
    echo ""
    ALL_NEW=0
    if [ $NEW_PAR -eq 2 ]; then
      ALL_NEW=1
    fi
    $THIS/distTree_inc_new.sh $INC $ALL_NEW
    if [ -e $INC/finished ]; then
      break
    fi
  done
fi
  

VARIANCE=`cat $INC/variance`


echo ""
echo ""
echo "Final optimization ..."
VER=`cat $INC/version`
echo "$VER  # Final optimization  `date`  `date +%s`" >> $INC/runlog  
cp $INC/tree $INC/hist/tree.$VER
gzip $INC/hist/tree.$VER
#
VER=$(( $VER + 1 ))
echo $VER > $INC/version


DELETE=""
if [ -e $INC/outlier-genogroup ]; then
  wc -l $INC/outlier-genogroup
  DELETE="-delete $INC/outlier-genogroup  -check_delete"
fi

# Time: O(n log^4(n))
# PAR
$THIS/makeDistTree  -threads 30  -data $INC/  -variance $VARIANCE  $DELETE \
  -optimize  -skip_len  -subgraph_iter_max 5 \
  -output_tree $INC/tree.new  -leaf_errors leaf_errors  > $INC/hist/makeDistTree-complete-inc.$VER
mv $INC/tree.new $INC/tree
# -reinsert  
#tail -n +5 leaf_errors.dm | sort -k 2 -g -r > leaf_errors.txt

if [ -e $INC/outlier-genogroup ]; then
  echo ""
  echo "Database: genogroup outliers ..."
  $INC/objects_in_tree.sh $INC/outlier-genogroup null
  mv $INC/outlier-genogroup $INC/hist/outlier-genogroup.$VER
fi


echo ""
echo "QC ..."
$INC/qc.sh go
echo ""
echo "Tree QC ..."
$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE  -qc  -noqual > $INC/hist/makeDistTree-qc.$VER
else
  VER=`cat $INC/version`
  VARIANCE=`cat $INC/variance`
fi 


$THIS/distTree_inc_tree1_quality.sh $INC



if [ -e $INC/phen ]; then
	DATE=`date +%Y%m%d`
	PHEN_LARGE=`cat $INC/phen_large`

	echo ""
	echo "Root and quality ..."
	$THIS/tree_quality_phen.sh $INC/tree "" $INC/phen $PHEN_LARGE 1 "" qual.raw > $INC/hist/tree_quality_phen.$VER 
	cat $INC/hist/tree_quality_phen.$VER 
	OLD_ROOT=`grep '^Old root: ' $INC/hist/tree_quality_phen.$VER | sed 's/^Old root: //1'`
	NEW_ROOT=`grep '^New root: ' $INC/hist/tree_quality_phen.$VER | sed 's/^New root: //1'`

	echo ""
	echo "Setting root and sorting ..."
  if [ ! "$NEW_ROOT" ]; then
    NEW_ROOT=$OLD_ROOT
  fi
  # -noqual must be absent to compute quality data after reroot()
	$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE  -reroot_at "$NEW_ROOT"  -output_tree tree.$DATE > /dev/null
	
	echo ""
	echo "Names ..."
	$THIS/tree2names.sh tree.$DATE $INC/phen $PHEN_LARGE > $INC/hist/tree2names.$VER
fi


echo ""
date

