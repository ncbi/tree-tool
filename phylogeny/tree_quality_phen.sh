#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Phenotypic quality of a distance tree, find root"
  echo "#1: distance tree"
  echo "#2: target list of objects (superset) | '' - all"
  echo "#3: phen/"
  echo "#4: large directories (0/1)"
  echo "#5: find root, save RAM (0/1)"
  echo "#6: output .qual-file or ''"
  echo "Time: 3.5 hours/172K genomes"
  exit 1
fi
TREE=$1
TARGET="$2"
PHEN=$3
LARGE=$4
FIND_ROOT=$5
QUAL="$6"


$THIS/../check_tmp.sh
TMP=`mktemp`
comment $TMP


DELETE=""
if [ "$TARGET" ]; then
  sort -cu $TARGET
  $THIS/tree2obj.sh $TREE > $TMP.cur
  $THIS/../setMinus $TMP.cur $TARGET > $TMP.del
  DELETE="-delete $TMP.del  -check_delete"
fi
$THIS/makeDistTree  -threads 15  -input_tree $TREE  $DELETE  -noqual  -output_feature_tree $TMP.feature_tree > $TMP.distTree

LARGE_PAR=""
if [ $LARGE == 1 ]; then
  LARGE_PAR="-large"
fi
FIND_ROOT_PARAM=""
if [ $FIND_ROOT == 1 ]; then
  FIND_ROOT_PARAM="-output_core $TMP.core  -save_mem"
fi
$THIS/makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree  -features $PHEN  $LARGE_PAR  -nominal_singleton_is_optional  $FIND_ROOT_PARAM  -qual $TMP.qual
  # -prefer_gain  

if [ "$QUAL" ]; then
  cp $TMP.qual $QUAL
fi


rm -f $TMP*
