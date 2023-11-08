#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Build a distance tree incrementally"
  echo "Update: #1/"
  echo "Output: leaf_errors.dm.gz, arc_existence.dm.gz"
  echo "        if #1/phen exists then: tree.<DATE>, disagreement_nodes[.txt], disagreement_objects, gain_nodes, qual, qual.raw"
  echo "Requires: large RAM, large running time"
  echo "#1: incremental distance tree directory"
  echo "#2: add new objects from #1/new/ (0/1)"
  echo "#3: list of good quality objects or ''"
  echo "#4: release directory where subdirectories are numbers, or ''. Valid if #1/phen exists"
  echo "Time: O(n log^4(n))"
  exit 1
fi
INC=$1
NEW_PAR=$2
GOOD="$3"
RELDIR="$4"


$THIS/../check_tmp.sh


QC=1


if [ $QC == 1 ]; then
  section "QC"
  $INC/qc.sh 0
fi


VARIANCE=`cat $INC/variance`

N=15
if [ -e $INC/threads ]; then
  N=`cat $INC/threads`
fi
THREADS="-threads $N"


if true; then   
if [ $NEW_PAR == 1 ]; then
  # Time: O(n log^4(n))
  while true; do
    if [ -e $INC/stop ]; then
      warning '*** STOPPED ***'
      exit 2
    fi
    
    if [ -e $INC/skip ]; then
      warning '*** SKIPPED ***'
      break
    fi
    
    VER=`cat $INC/version`
    echo "$VER  `date`  `date +%s`" >> $INC/runlog  
    pwd
    echo ""
    $THIS/distTree_inc_new.sh $INC 
    if [ -e $INC/finished ]; then
      break
    fi
  done
fi
  

super_section "Final optimization"

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
$THIS/makeDistTree  $THREADS  -data $INC/  -variance $VARIANCE  $DELETE \
  -optimize  -skip_len  -subgraph_iter_max 2 \
  -output_tree $INC/tree.new  -leaf_errors leaf_errors  -arc_existence arc_existence  > $INC/hist/makeDistTree-final.$VER
mv $INC/tree.new $INC/tree
# -reinsert  
#tail -n +5 leaf_errors.dm | sort -k 2 -g -r > leaf_errors.txt

section "arc_existence probabilities"
tail -n +6 arc_existence.dm | awk '{print $3};'| $THIS/../dm/count | grep -w "count\|sum"

gzip leaf_errors.dm
gzip arc_existence.dm

if [ -e $INC/outlier-genogroup ]; then
  section "Database: genogroup outliers"
  $INC/objects_in_tree.sh $INC/outlier-genogroup null
  mv $INC/outlier-genogroup $INC/hist/outlier-genogroup.$VER
fi


if [ $QC == 1 ]; then
  section "QC"
  $INC/qc.sh 0
fi
section "Tree QC"
$THIS/makeDistTree  $THREADS  -data $INC/  -variance $VARIANCE  -qc  -noqual > $INC/hist/makeDistTree-qc.$VER
else
  VER=`cat $INC/version`
fi 


$THIS/distTree_inc_tree1_quality.sh $INC  


if [ -e $INC/phen ]; then
	DATE=`date +%Y%m%d`

	LARGE=0
	if [ -e $INC/large ]; then
	  LARGE=1
	fi

	super_section "Root and quality of whole tree"
	$THIS/tree_quality_phen.sh $INC/tree "" $INC/phen $LARGE 1 qual.raw > $INC/hist/tree_quality_phen.$VER 
	cat $INC/hist/tree_quality_phen.$VER 
	OLD_ROOT=`grep '^Old root: ' $INC/hist/tree_quality_phen.$VER | sed 's/^Old root: //1'`
	NEW_ROOT=`grep '^New root: ' $INC/hist/tree_quality_phen.$VER | sed 's/^New root: //1'`

	section "Setting root and sorting"
  if [ ! $NEW_ROOT ]; then
    NEW_ROOT=$OLD_ROOT
  fi
  # -noqual must be absent to compute quality data after reroot()
	$THIS/makeDistTree  $THREADS  -data $INC/  -variance $VARIANCE  -reroot_at "$NEW_ROOT"  -output_tree tree.$DATE > /dev/null
	
	super_section "Names"
	$THIS/tree2names.sh tree.$DATE $INC/phen $LARGE > $INC/hist/tree2names.$VER

  if [ "$GOOD" ]; then
  	super_section "Quality of good quality tree"
  	$THIS/tree_quality_phen.sh $INC/tree $GOOD $INC/phen $LARGE 1 "" > $INC/hist/tree_quality_phen-good.$VER 
  	cat $INC/hist/tree_quality_phen-good.$VER 
  fi


  if [ -n "$RELDIR" ]; then
    super_section "Release"
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
    mv disagreement_objects disagreement_nodes.txt disagreement_nodes gain_nodes qual tree.$DATE qual.raw leaf_errors.dm.gz arc_existence.dm.gz $RELDIR/$RELNUM/
    
    rm -f $RELDIR/latest
    ln -s $PWD/$RELDIR/$RELNUM $RELDIR/latest
    
    rm -f $INC/tree.released
    ln -s $PWD/$RELDIR/$RELNUM/tree.$DATE $INC/tree.released
  fi
fi


success

