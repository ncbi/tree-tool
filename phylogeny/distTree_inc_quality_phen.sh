#!/bin/bash --noprofile
THIS=$( dirname $0 )
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
INC=$1
START=$2 
TARGET=$3
VER_INC=$4
SEED=$5


[ $START -ge 1 ]       ||  exit 1
[ $START -le $TARGET ] ||  exit 1
[ $VER_INC -ge 1 ]     ||   exit 1
[ $SEED -gt 0 ]        ||  exit 1


TMP=$( mktemp )
comment $TMP


$THIS/tree2obj.sh $INC/hist/tree.$START > $TMP.start
$THIS/tree2obj.sh $INC/hist/tree.1 > $TMP.init
$THIS/../setMinus $TMP.start $TMP.init > $TMP.start-new
$THIS/tree2obj.sh $INC/tree > $TMP.final
$THIS/../setIntersect.sh $TMP.start-new $TMP.final 0 > $TMP.start-good

N=$( < $TMP.start-good  wc -l )
echo "Tree size: $N"
S=$(( N / 2 ))
S_max=300  # PAR
if [ $S -gt $S_max ]; then
  S=$S_max
fi
echo "Sample size: $S"

$THIS/../setRandOrd $TMP.start-good  -seed $SEED  -sigpipe | head -$S > $TMP.test
if [ ! -s $TMP.test ]; then
	cp $TMP.init $TMP.test
  cp $INC/hist/tree.$START $TMP.tree
	$THIS/makeDistTree  -input_tree $TMP.tree  -output_feature_tree $TMP.feature_tree1 | grep "# Discernible leaves:"
else
	section "Complete optimization"
	cat $TMP.init >> $TMP.test
	$THIS/../sort.sh $TMP.test
	$THIS/../list2pairs $TMP.test > $TMP.pairs
	$THIS/distTree_inc_request2dissim.sh $INC $TMP.pairs $TMP.dissim
	$THIS/../dm/conversion/pairs2dm $TMP.dissim 1 cons 6 -distance > $TMP.dm
	$THIS/makeDistTree  -data $TMP  -dissim_attr cons  -optimize  -output_tree $TMP.tree  -output_feature_tree $TMP.feature_tree1 | grep "# Discernible leaves:"
fi

$THIS/makeFeatureTree  -input_tree $TMP.feature_tree1  -features $INC/phen  -nominal_singleton_is_optional  -qual $TMP.qual  


mkdir $TMP.trees
i=$START
n=0
while [ $i -le $TARGET ]
do
	$THIS/distTree_inc_quality_phen_step.sh $INC $TMP.test $i $TMP
	i=$(( i + VER_INC ))
	n=$(( n + 1 ))
done

$THIS/distTree_inc_quality_phen_step.sh $INC $TMP.test opt $TMP


rm -fr $TMP*
