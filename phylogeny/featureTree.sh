#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 4 ]; then
  echo "Output: #3[.maxParsimony].{tree,core} - feature tree"
  echo "#1: list of genomes"
  echo "#2: directory with features of genomes"
  echo "#3: use time (1/0)"
  echo "#4: output prefix"
  exit 1
fi
OBJ=$1
FEAT=$2
USE_TIME=$3
OUT=$4


TMP=`mktemp`
comment $TMP
#set -x


QC="" # -qc


$THIS/../dissim/feature2dissim $QC $OBJ $FEAT > $TMP.dm

super_section "Building distance tree"
$THIS/makeDistTree  $QC  -data $TMP  -dissim_attr dissim  -variance linExp  -optimize  -subgraph_iter_max 20  -output_feature_tree ${TMP}-init.tree 

super_section "Building feature tree without time"
$THIS/makeFeatureTree  $QC  -input_tree ${TMP}-init.tree  -features $FEAT  -optim_iter_max 100  -prefer_gain  -output_core $OUT.maxParsimony.core  -output_tree $OUT.maxParsimony.tree  

if [ $USE_TIME == 1 ]; then
  # MLE
  super_section "Building feature tree with time"
  $THIS/makeFeatureTree  $QC  -input_tree $OUT.maxParsimony.tree  -features $FEAT  -input_core $OUT.maxParsimony.core  -optim_iter_max 100  -use_time  -output_tree $OUT.tree  -output_core $OUT.core 
fi


rm $TMP*
