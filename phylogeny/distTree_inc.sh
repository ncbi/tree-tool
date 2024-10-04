#!/bin/bash --noprofile
THIS=$( dirname $0 )
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
NEW=$2
GOOD="$3"
RELDIR="$4"

QC=1


$THIS/../check_tmp.sh

if [ $GOOD ]; then
  $THIS/../check_file.sh $GOOD 1
fi

if [ $QC == 1 ]; then
  section "QC"
  $INC/qc.sh 0
fi



if true; then   
  if [ $NEW == 1 ]; then
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
      
      VER=$( cat $INC/version )
      echo "$VER  $( date )  $( date +%s )" >> $INC/runlog  
      pwd
      echo ""
      $THIS/distTree_inc_new.sh $INC 
      if [ -e $INC/finished ]; then
        break
      fi
    done
  fi
    

  super_section "Final optimization"

  VER=$( cat $INC/version )
  echo "$VER  # Final optimization  $( date)   $( date +%s )" >> $INC/runlog  
  cp $INC/tree $INC/hist/tree.$VER
  if [ $VER -gt 1 ]; then
    gzip $INC/hist/tree.$VER
  fi
  #
  VER=$(( $VER + 1 ))
  echo $VER > $INC/version

  VARIANCE=$( cat $INC/variance )

  # THREADS
  N=15
  if [ -e $INC/threads ]; then
    N=$( cat $INC/threads )
  fi
  THREADS="-threads $N"

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
  VER=$( cat $INC/version )
fi 


$THIS/distTree_inc_tree1_quality.sh $INC  


if [ -e $INC/phen ]; then
	LARGE=0
	if [ -e $INC/large ]; then
	  LARGE=1
	fi

  DATE=$( date +%Y%m%d )
  OUT_TREE=tree.$DATE
  
  $THIS/distTree_release.sh $INC/tree $INC/phen $LARGE 1 $OUT_TREE $RELDIR 
  
  if [ "$GOOD" ]; then
  	super_section "Quality of good quality tree"
  	$THIS/tree_quality_phen.sh $INC/tree $GOOD $INC/phen $LARGE 1 "" > $INC/hist/tree_quality_phen-good.$VER 
  	cat $INC/hist/tree_quality_phen-good.$VER 
  fi

  if [ -n "$RELDIR" ]; then
    LATEST=$( realpath $PWD/$RELDIR/latest )
    mv leaf_errors.dm.gz arc_existence.dm.gz $LATEST
    rm -f $INC/tree.released
    ln -s $LATEST/$OUT_TREE $INC/tree.released
  fi
fi


success

