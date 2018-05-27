#!/bin/csh -f

if ($# != 2) then
  echo "Input: sm1K-univ-calibrate.dm"
  echo "#1: power"
  echo "#2: phen/"
  exit 1
endif


echo ""
echo "power = $1"

calibrateDissims sm1K-univ-calibrate $1  -output_dissim sm1K-univ.pairs  > hmm-univ.stat
if ($?) exit 1

tail -n +5 sm1K-univ.pairs.dm | sed 's/-/ /1' > sm1K-univ.pairs
if ($?) exit 1
rm sm1K-univ.pairs.dm

pairs2attr2 sm1K-univ.pairs 1 cons 6 -distance > sm1K-univ.dm
if ($?) exit 1
rm sm1K-univ.pairs

distTree.sh sm1K-univ cons
if ($?) exit 1
rm sm1K-univ.dm

tree_quality_phen.sh sm1K-univ.tree $2
if ($?) exit 1


