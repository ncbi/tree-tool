#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 2 ]; then
  echo "Output: #1.{tree,core} - feature tree"
  echo "#1: #1.list - list of genomes"
  echo "#2: directory with features of genomes"
  exit 1
fi



$THIS/feature2dissim $1.list $2 > $1.dm

$THIS/distTree.sh $1 cons 
rm $1.dm
rm $1.makeDistTree

# Gain/loss tree
$THIS/makeDistTree  -input_tree $1.tree  -noqual  -output_feature_tree $1-init.tree  
rm $1.tree 

echo ""
echo ""
$THIS/makeFeatureTree  -input_tree $1-init.tree  -features $2  -output_tree $1-maxParsimony.tree  -output_core $1-maxParsimony.core  -optim_iter_max 100
rm $1-init.tree

# MLE
echo ""
echo ""
$THIS/makeFeatureTree  -input_tree $1-maxParsimony.tree  -features $2  -input_core $1-maxParsimony.core  -output_tree $1.tree  -output_core $1.core  -optim_iter_max 100  -use_time
rm $1-maxParsimony.tree 
rm $1-maxParsimony.core


# To distance tree for loading into Phyl
# mv genomes299.tree genomes299.featureTree
# tail -n +3 genomes299.featureTree | sed 's/: t=/: len=/1' | grep -v '^ *g' | sed 's/^\( *\)s/\1/1' | sed 's/^1: len=0.000000e+00/1: len=nan/1' | sed 's/^\( *\)\([1-9]\)/\10x\2/1' > genomes299.tree
