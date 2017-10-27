#!/bin/csh -f

if ($# != 2) then
  echo "Output: #1.dm - distances, #1.distTree, #1-maxParsimony.{tree,core}, #1.{tree,core} - feature tree"
  echo "#1: #1.list - list of genomes"
  echo "#2: directory with features of genomes"
  exit 1
endif



feature2dist $1.list $2 > $1.dm
if ($?) exit 1

# Add comments to $1.dm

distTree.sh $1 cons 
if ($?) exit 1 

# Gain/loss tree
makeDistTree  -input_tree $1.tree  -output_feature_tree $1.featureTree
if ($?) exit 1 
mv $1.tree $1.distTree  

# Make all times = 0
cat $1.featureTree | sed 's/: t=\([^ ]\+\) /: t=nan /1' > $1-init.tree
if ($?) exit 1 
rm $1.featureTree

makeFeatureTree  -input_tree $1-init.tree  -genes $2  -output_tree $1-maxParsimony.tree  -output_core $1-maxParsimony.core  -optim_iter_max 100
if ($?) exit 1 
rm $1-init.tree

# MLE
makeFeatureTree  -input_tree $1-maxParsimony.tree  -genes $2  -input_core $1-maxParsimony.core  -output_tree $1.tree  -output_core $1.core  -optim_iter_max 100  -use_time
if ($?) exit 1 



# To distance tree for loading into Phyl
# mv genomes299.tree genomes299.featureTree
# tail -n +3 genomes299.featureTree | sed 's/: t=/: len=/1' | grep -v '^ *g' | sed 's/^\( *\)s/\1/1' | sed 's/^1: len=0.000000e+00/1: len=nan/1' | sed 's/^\( *\)\([1-9]\)/\10x\2/1' > genomes299.tree
