#!/bin/csh -f

if ($# != 2) then
  echo "Input: sm1K.combo, dissim_scale"
  echo "#1: coeff (> 0)"
  echo "#2: phen/"
  exit 1
endif


combine_dissims sm1K.combo dissim_scale $1 > sm1K.dissim
if ($?) exit 1

pairs2attr2 sm1K.dissim 1 cons 6 -distance > sm1K.dm
if ($?) exit 1
rm sm1K.dissim

echo ""
distTree.sh sm1K cons
if ($?) exit 1
rm sm1K.dm

echo ""
echo "Tree quality ..."
tree_quality_phen.sh sm1K.tree $2
if ($?) exit 1

