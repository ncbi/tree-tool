#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Build a distance tree incrementally"
  echo "Update: #1/"
  echo "Output: leaf_errors.dm"
  echo "        if #1/phen exists then: tree.<DATE>, disagreement_nodes[.txt], disagreement_objects, gain_nodes, qual, qual.raw"
  echo "Requires: large RAM, large running time"
  echo "#1: incremental distance tree directory"
  echo "#2: add new objects from #1/new/ (0/1)"
  echo "#3: release directory where subdirectories are numbers, or ''. Valid if #1/phen exists"
  echo "Time: O(n log^4(n))"
  exit 1
fi
INC=$1
NEW_PAR=$2
RELDIR="$3"


QC=1


if [ $QC == 1 ]; then
  echo "QC ..."
  $INC/qc.sh go
  echo ""
fi


if true; then   
if [ $NEW_PAR == 1 ]; then
  # Time: O(n log^4(n))
  while true; do
    if [ -e $INC/stop ]; then
      echo ""
      echo -e "$YELLOW"'*** STOPPED ***'"$NOCOLOR"
      exit 2
    fi
    
    if [ -e $INC/skip ]; then
      echo ""
      echo -e "$YELLOW"'*** SKIPPED ***'"$NOCOLOR"
      break
    fi
    
    NEW=`ls $INC/new/ | wc -l`
    VER=`cat $INC/version`
    echo "$VER  # New: $NEW  `date`  `date +%s`" >> $INC/runlog  
    echo ""
    echo ""
    $THIS/distTree_inc_new.sh $INC 
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
if [ $VER -gt 1 ]; then
  gzip $INC/hist/tree.$VER
fi
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


if [ $QC == 1 ]; then
  echo ""
  echo "QC ..."
  $INC/qc.sh go
fi
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
	$THIS/tree_quality_phen.sh $INC/tree "" $INC/phen $PHEN_LARGE 1 qual.raw > $INC/hist/tree_quality_phen.$VER 
	cat $INC/hist/tree_quality_phen.$VER 
	OLD_ROOT=`grep '^Old root: ' $INC/hist/tree_quality_phen.$VER | sed 's/^Old root: //1'`
	NEW_ROOT=`grep '^New root: ' $INC/hist/tree_quality_phen.$VER | sed 's/^New root: //1'`

	echo ""
	echo "Setting root and sorting ..."
  if [ ! $NEW_ROOT ]; then
    NEW_ROOT=$OLD_ROOT
  fi
  # -noqual must be absent to compute quality data after reroot()
	$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE  -reroot_at "$NEW_ROOT"  -output_tree tree.$DATE > /dev/null
	
	echo ""
	echo "Names ..."
	$THIS/tree2names.sh tree.$DATE $INC/phen $PHEN_LARGE > $INC/hist/tree2names.$VER


  if [ -n "$RELDIR" ]; then
    echo ""
    echo ""
    if [ ! -e $RELDIR ]; then
      error "$RELDIR does not exist"
    fi
    set +o errexit
    RELNUM=`ls $RELDIR | grep -v '[^0-9]' | sort -n -r | head -1`
    set -o errexit
    if [ -z "$RELNUM" ]; then
      RELNUM=0
    fi
    RELNUM=$(( $RELNUM + 1 ))
    echo "Release: $RELNUM"

    mkdir $RELDIR/$RELNUM
    mv disagreement_objects disagreement_nodes.txt disagreement_nodes gain_nodes qual tree.$DATE qual.raw leaf_errors.dm $RELDIR/$RELNUM/
    rm -f $RELDIR/latest
    ln -s $PWD/$RELDIR/$RELNUM $RELDIR/latest
    
    rm -f $INC/tree.released
    ln -s $PWD/$RELDIR/$RELNUM/tree.$DATE $INC/tree.released
  fi
fi


echo ""
date

