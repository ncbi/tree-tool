#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Phenotypic quality comparison of 2 distance trees"
  echo "Output: qual.comp, gain1, gain2, loss1, loss2"
  echo "#1: distance tree 1"
  echo "#2: distance tree 2"
  echo "#3: phen/"
  echo "#4: phen/ is a large directory (0/1)"
  echo "#5: list of objects to intersect with (sorted) | ''"
  exit 1
fi
T1=$1
T2=$2
PHEN=$3
LARGE=$4
TARGET="$5"


$THIS/../check_tmp.sh

TMP=$( mktemp )
comment $TMP 


echo "$T1 vs. $T2"


$THIS/tree2obj.sh $T1 > $TMP.list1
$THIS/tree2obj.sh $T2 > $TMP.list2

# QC
uniq $TMP.list1 > $TMP.list1-uniq
uniq $TMP.list2 > $TMP.list2-uniq
diff $TMP.list1 $TMP.list1-uniq
diff $TMP.list2 $TMP.list2-uniq

if [ "$TARGET" ]; then
  wc -l $TARGET
  sort -cu $TARGET
  $THIS/../setIntersect.sh $TMP.list1 $TARGET 0 > $TMP.list1-good 
  $THIS/../setIntersect.sh $TMP.list2 $TARGET 0 > $TMP.list2-good 
else
  cp $TMP.list1 $TMP.list1-good
  cp $TMP.list2 $TMP.list2-good
fi
wc -l $TMP.list1-good
wc -l $TMP.list2-good

N=$( < $TMP.list1-good  wc -l )
if [ $N -le 2 ]; then
  exit 0
fi

N=$( < $TMP.list2-good  wc -l )
if [ $N -le 2 ]; then
  exit 0
fi

N=$( $THIS/../setIntersect.sh $TMP.list1-good $TMP.list2-good 0 | wc -l )
if [ $N -le 2 ]; then
  echo "Intersection = $N"
  exit 0
fi

$THIS/../setMinus $TMP.list1 $TMP.list2-good > $TMP.list1-del
$THIS/../setMinus $TMP.list2 $TMP.list1-good > $TMP.list2-del

$THIS/makeDistTree  -threads 15  -input_tree $T1  -delete $TMP.list1-del  -check_delete  -noqual  -output_tree $TMP.tree1  -output_feature_tree $TMP.feature_tree1 > $TMP.makeDistTree1
$THIS/makeDistTree  -threads 15  -input_tree $T2  -delete $TMP.list2-del  -check_delete  -noqual  -output_tree $TMP.tree2  -output_feature_tree $TMP.feature_tree2 > $TMP.makeDistTree2


LARGE_PAR=""
if [ $LARGE == 1 ]; then
  LARGE_PAR="-large"
fi

section "$T1"
$THIS/makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree1  -features $PHEN  $LARGE_PAR  -nominal_singleton_is_optional  -output_core $TMP.core1  -qual $TMP.qual1  -gain_nodes gain1  -loss_nodes loss1  -disagreement_nodes $TMP.disagreement_nodes1
  # -prefer_gain  
section "$T2"
$THIS/makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree2  -features $PHEN  $LARGE_PAR  -nominal_singleton_is_optional  -output_core $TMP.core2  -qual $TMP.qual2  -gain_nodes gain2  -loss_nodes loss2  -disagreement_nodes $TMP.disagreement_nodes2
  # -prefer_gain  
  
# .qual:
# name gains loss criterion objects optional_objects
# 1    2     3    4         5       6

# QC
cut -f 1 $TMP.qual1 > $TMP.qual1.list
cut -f 1 $TMP.qual2 > $TMP.qual2.list
diff $TMP.qual1.list $TMP.qual2.list

echo -e "#Feature\ttau1_less_tau2" > qual.comp
join -1 1 -2 1 -t$'\t' $TMP.qual1 $TMP.qual2 | awk -F '\t' '{OFS="\t"; D=($2 + $3) - ($7 + $8); if (D != 0) print $1, D};' >> qual.comp

    DIFF=$( tail -n +2 qual.comp | cut -f 2 |               count | grep -w '^sum' | cut -f 2 )
DIFF_NEG=$( tail -n +2 qual.comp | cut -f 2 | grep    '^-'| count | grep -w '^sum' | cut -f 2 )
DIFF_POS=$( tail -n +2 qual.comp | cut -f 2 | grep -v '^-'| count | grep -w '^sum' | cut -f 2 )
echo "Difference: $DIFF = $DIFF_POS $DIFF_NEG"


rm -f $TMP*
