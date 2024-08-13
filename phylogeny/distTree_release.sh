#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Release a distance tree"
  echo "Output: #5, disagreement_nodes[.txt], disagreement_objects, gain_nodes, qual, qual.raw, tree_quality_phen.out, tree2names.out"
  echo "Requires: large RAM, large running time"
  echo "#1: tree"
  echo "#2: phen/"
  echo "#3: #2 is large (0/1)"
  echo "#4: reroot (0/1)"
  echo "#5: output tree file name"
  echo "#6: release directory where subdirectories are numbers | ''"
  exit 1
fi
IN_TREE=$1
PHEN=$2
LARGE=$3
REROOT=$4
OUT_TREE=$5
RELDIR=$6


super_section "Root and quality of whole tree"
$THIS/tree_quality_phen.sh $IN_TREE "" $PHEN $LARGE 1 qual.raw > tree_quality_phen.out
cat tree_quality_phen.out

# $OUT_TREE
if [ $REROOT == 1 ]; then
  section "Setting root and sorting"
  OLD_ROOT=$( grep '^Old root: ' tree_quality_phen.out | sed 's/^Old root: //1' )
  NEW_ROOT=$( grep '^New root: ' tree_quality_phen.out | sed 's/^New root: //1' )
  if [ ! $NEW_ROOT ]; then
    NEW_ROOT=$OLD_ROOT
  fi
  # -noqual must be absent to compute quality data after reroot()
  #$THIS/makeDistTree  -data $INC/  -variance $VARIANCE  -reroot_at "$NEW_ROOT"  -output_tree $OUT_TREE > /dev/null
  $THIS/makeDistTree  -input_tree $IN_TREE  -reroot_at "$NEW_ROOT"  -output_tree $OUT_TREE  > /dev/null
else
  cp $IN_TREE $OUT_TREE
fi

super_section "Names"
$THIS/tree2names.sh $OUT_TREE $PHEN $LARGE > tree2names.out


if [ -n "$RELDIR" ]; then
  super_section "Release"

  if [ ! -e $RELDIR ]; then
    error "$RELDIR does not exist"
  fi

  # RELNUM  
  set +o errexit
  RELNUM=`ls $RELDIR | grep -v '[^0-9]' | sort -n -r | head -1`
  set -o errexit
  if [ -z "$RELNUM" ]; then
    RELNUM=0
  fi
  RELNUM=$(( $RELNUM + 1 ))
  echo "Release: $RELNUM"

  mkdir $RELDIR/$RELNUM
  mv disagreement_objects disagreement_nodes.txt disagreement_nodes gain_nodes qual $OUT_TREE qual.raw tree_quality_phen.out tree2names.out $RELDIR/$RELNUM/
  
  rm -f $RELDIR/latest
  ln -s $PWD/$RELDIR/$RELNUM $RELDIR/latest
fi




