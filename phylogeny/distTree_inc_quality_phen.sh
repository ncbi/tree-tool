#!/bin/bash
source bash_common.sh
if [ $# -ne 5 ]; then
  echo "Phenotypic quality of an incremental tree"
  echo "Requires: no added outliers"
  echo "#1: incremental tree"
  echo "#2: start tree version"
  echo "#3: target tree version (>= #2)"
  echo "#4: version increment (>= 1)"
  echo "#5: seed (>= 1)"
  exit 1
fi


if ($2 < 1)   exit 1
if ($2 > $3)  exit 1
if ($4 < 1)   exit 1
if ($5 <= 0)  exit 1


TMP=`mktemp`
echo $TMP


tree2obj.sh $1/hist/tree.$2 > $TMP.start
tree2obj.sh $1/hist/tree.1 > $TMP.init
setMinus $TMP.start $TMP.init > $TMP.start-new
tree2obj.sh $1/tree > $TMP.final
setIntersect.sh $TMP.start-new $TMP.final 0 > $TMP.start-good

N=`cat $TMP.start-good | wc -l`
echo "Tree size: $N"
S=$(( $N / 2 ))
S_max=300  # PAR
if [ $S > $S_max ]; then
  S=$S_max
fi
echo "Sample size: $S"

setRandOrd $TMP.start-good  -seed $5  -sigpipe | head -$S > $TMP.test
if [ ! -s $TMP.test ]; then
	cp $TMP.init $TMP.test
  cp $1/hist/tree.$2 $TMP.tree
	makeDistTree  -input_tree $TMP.tree  -output_feature_tree $TMP.feature_tree1 | grep "# Discernible leaves:"
else
	echo ""
	echo "Complete optimization  ..."
	cat $TMP.init >> $TMP.test
	sort.sh $TMP.test
	list2pairs $TMP.test > $TMP.pairs
	distTree_inc_request2dissim.sh $1 $TMP.pairs $TMP.dissim
	pairs2attr2 $TMP.dissim 1 cons 6 -distance > $TMP.dm
	makeDistTree  -data $TMP  -dissim cons  -optimize  -output_tree $TMP.tree  -output_feature_tree $TMP.feature_tree1 | grep "# Discernible leaves:"
fi

makeFeatureTree  -input_tree $TMP.feature_tree1  -features $1/phen  -nominal_singleton_is_optional  -output_core $TMP.core  -qual $TMP.qual  


mkdir $TMP.trees
i=$2
n=0
while [ $i <= $3 ]
do
	distTree_inc_quality_phen_step.sh $1 $TMP.test $i $TMP
	i=$(( $i + $4 ))
	n=$(( $n + 1 ))
done

distTree_inc_quality_phen_step.sh $1 $TMP.test opt $TMP

rm -fr $TMP*
