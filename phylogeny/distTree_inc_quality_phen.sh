#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
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


$THIS/tree2obj.sh $1/hist/tree.$2 > $TMP.start
$THIS/tree2obj.sh $1/hist/tree.1 > $TMP.init
$THIS/../setMinus $TMP.start $TMP.init > $TMP.start-new
$THIS/tree2obj.sh $1/tree > $TMP.final
$THIS/../setIntersect.sh $TMP.start-new $TMP.final 0 > $TMP.start-good

N=`cat $TMP.start-good | wc -l`
echo "Tree size: $N"
S=$(( $N / 2 ))
S_max=300  # PAR
if [ $S > $S_max ]; then
  S=$S_max
fi
echo "Sample size: $S"

$THIS/../setRandOrd $TMP.start-good  -seed $5  -sigpipe | head -$S > $TMP.test
if [ ! -s $TMP.test ]; then
	cp $TMP.init $TMP.test
  cp $1/hist/tree.$2 $TMP.tree
	$THIS/makeDistTree  -input_tree $TMP.tree  -output_feature_tree $TMP.feature_tree1 | grep "# Discernible leaves:"
else
	echo ""
	echo "Complete optimization  ..."
	cat $TMP.init >> $TMP.test
	sort.sh $TMP.test
	$THIS/../list2pairs $TMP.test > $TMP.pairs
	$THIS/distTree_inc_request2dissim.sh $1 $TMP.pairs $TMP.dissim
	$THIS/../dm/pairs2attr2 $TMP.dissim 1 cons 6 -distance > $TMP.dm
	$THIS/makeDistTree  -data $TMP  -dissim cons  -optimize  -output_tree $TMP.tree  -output_feature_tree $TMP.feature_tree1 | grep "# Discernible leaves:"
fi

$THIS/makeFeatureTree  -input_tree $TMP.feature_tree1  -features $1/phen  -nominal_singleton_is_optional  -qual $TMP.qual  


mkdir $TMP.trees
i=$2
n=0
while [ $i <= $3 ]
do
	$THIS/distTree_inc_quality_phen_step.sh $1 $TMP.test $i $TMP
	i=$(( $i + $4 ))
	n=$(( $n + 1 ))
done

$THIS/distTree_inc_quality_phen_step.sh $1 $TMP.test opt $TMP


rm -fr $TMP*
