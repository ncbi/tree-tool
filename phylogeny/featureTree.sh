#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 2 ]; then
  echo "Output: #1.{tree,core} - feature tree"
  echo "#1: #1.list - list of genomes"
  echo "#2: directory with features of genomes"
  exit 1
fi
OBJ=$1
FEAT=$2


TMP=`mktemp`


QC="" # -qc


$THIS/../dissim/feature2dissim $QC $OBJ.list $FEAT > $TMP.dm

echo ""
echo ""
$THIS/makeDistTree  $QC  -data $TMP  -dissim_attr dissim  -variance linExp  -optimize  -output_feature_tree ${TMP}-init.tree 

echo ""
echo ""
$THIS/makeFeatureTree  $QC  -input_tree ${TMP}-init.tree  -features $FEAT  -optim_iter_max 100  -output_core ${TMP}-maxParsimony.core  -prefer_gain  -output_tree ${TMP}-maxParsimony.tree  

# MLE
echo ""
echo ""
$THIS/makeFeatureTree  $QC  -input_tree ${TMP}-maxParsimony.tree  -features $FEAT  -input_core ${TMP}-maxParsimony.core  -optim_iter_max 100  -use_time  -output_tree $OBJ.tree  -output_core $OBJ.core 


# Creating a distance tree for loading into Phyl:
# mv genomes299.tree genomes299.featureTree
# tail -n +3 genomes299.featureTree | sed 's/: t=/: len=/1' | grep -v '^ *g' | sed 's/^\( *\)s/\1/1' | sed 's/^1: len=0.000000e+00/1: len=nan/1' | sed 's/^\( *\)\([1-9]\)/\10x\2/1' > genomes299.tree


rm $TMP*
