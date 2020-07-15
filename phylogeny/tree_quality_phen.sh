  #!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Phenotypic quality of a distance tree, find root"
  echo "#1: distance tree"
  echo "#2: target list of objects | '' - all"
  echo "#3: phen/"
  echo "#4: phen_large (0/1)"
  echo "#5: find root and save RAM (0/1)"
  echo "#6: output .qual-file or ''"
  echo "Time: 3.5 hours/172K genomes"
  exit 1
fi
TREE=$1
TARGET="$2"
PHEN=$3
PHEN_LARGE=$4
FIND_ROOT=$5
QUAL="$6"


TMP=`mktemp`
echo $TMP


DELETE=""
if [ $TARGET ]; then
  $THIS/tree2obj.sh $TREE > $TMP.cur
  $THIS/../setMinus $TMP.cur $TARGET > $TMP.del
  DELETE="-delete $TMP.del  -check_delete"
fi
$THIS/makeDistTree  -threads 15  -input_tree $TREE  $DELETE  -noqual  -output_feature_tree $TMP.feature_tree > $TMP.distTree

LARGE=""
if [ $PHEN_LARGE == 1 ]; then
  LARGE="-large"
fi
FIND_ROOT_PARAM=""
if [ $FIND_ROOT == 1 ]; then
  FIND_ROOT_PARAM="-output_core $TMP.core  -save_mem"
fi
$THIS/makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree  -features $PHEN  $LARGE  -prefer_gain  -nominal_singleton_is_optional  $FIND_ROOT_PARAM  -qual $TMP.qual

if [ $QUAL ]; then
  cp $TMP.qual $QUAL
fi


rm -f $TMP*
