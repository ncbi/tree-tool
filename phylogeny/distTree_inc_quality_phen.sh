#!/bin/csh -f

if ($# != 5) then
  echo "Phenotypic quality of an incremental tree"
  echo "Requires: no added outliers"
  echo "#1: incremental tree"
  echo "#2: start tree version"
  echo "#3: target tree version (>= #2)"
  echo "#4: version increment (>= 1)"
  echo "#5: seed (>= 1)"
  exit 1
endif


if ($2 < 1)   exit 1
if ($2 > $3)  exit 1
if ($4 < 1)   exit 1
if ($5 <= 0)  exit 1


set tmp = `mktemp`
if ($?) exit 1
echo $tmp


tree2obj.sh $1/hist/tree.$2 > $tmp.start
if ($?) exit 1
tree2obj.sh $1/hist/tree.1 > $tmp.init
if ($?) exit 1
setMinus $tmp.start $tmp.init > $tmp.start-new
if ($?) exit 1
tree2obj.sh $1/tree > $tmp.final
if ($?) exit 1
setIntersect.sh $tmp.start-new $tmp.final 0 > $tmp.start-good
if ($?) exit 1

set N = `wc -l $tmp.start-good`
echo "Tree size: $N[1]"
@ S = $N[1] / 2
set S_max = 300  # PAR
if ($S > $S_max)  set S = $S_max
echo "Sample size: $S"

setRandOrd $tmp.start-good  -seed $5  -sigpipe | head -$S > $tmp.test
if (-z $tmp.test) then
	cp $tmp.init $tmp.test
	if ($?) exit 1
  cp $1/hist/tree.$2 $tmp.tree
	if ($?) exit 1
	makeDistTree  -input_tree $tmp.tree  -output_feature_tree $tmp.feature_tree1 | grep "# Discernible leaves:"
	if ($?) exit 1
else
	echo ""
	echo "Complete optimization  ..."
	cat $tmp.init >> $tmp.test
	sort.sh $tmp.test
	
	list2pairs $tmp.test > $tmp.pairs
	if ($?) exit 1
	
	distTree_inc_request2dissim.sh $1 $tmp.pairs $tmp.dissim
	if ($?) exit 1
	
	pairs2attr2 $tmp.dissim 1 cons 6 -distance > $tmp.dm
	if ($?) exit 1
	
	makeDistTree  -data $tmp  -dissim cons  -optimize  -output_tree $tmp.tree  -output_feature_tree $tmp.feature_tree1 | grep "# Discernible leaves:"
	if ($?) exit 1
endif

makeFeatureTree  -input_tree $tmp.feature_tree1  -features $1/phen  -output_core $tmp.core  -qual $tmp.qual  
if ($?) exit 1


mkdir $tmp.trees
if ($?) exit 1
set i = $2
set n = 0
while ($i <= $3)
	distTree_inc_quality_phen_step.sh $1 $tmp.test $i $tmp
	if ($?) exit 1
	@ i = $i + $4
	@ n = $n + 1
end

distTree_inc_quality_phen_step.sh $1 $tmp.test opt $tmp
if ($?) exit 1

rm -fr $tmp.*
