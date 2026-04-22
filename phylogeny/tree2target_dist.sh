#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
META="Metadata.tsv"
if [ $# -ne 4 ]; then
  echo "Print: <object> <closest target> <distance>"
  echo "Input: $META"
  echo "#1: tree"
  echo "#2: list of target objects"
  echo "#3: target column name"
  echo "#4: report"
  exit 1
fi
TREE=$1
TARGET=$2
TARGET_NAME="$3"
OUT=$4


$THIS/../check_file.sh $META 1


TMP=$( mktemp )
comment $TMP


$THIS/tree2obj.sh $TREE > $TMP.tree
sort -u $TARGET > $TMP.target_raw
$THIS/../setIntersect.sh $TMP.tree $TMP.target_raw 0 > $TMP.target
$THIS/../trav $TMP.target "sed 's/$/\t%f/1' $TMP.tree" | awk '$1 != $2' > $TMP.pair

section "statDistTree"
$THIS/statDistTree $TREE  -dist_request $TMP.pair  -dist_pairs $TMP.dist  -qc 
echo -e "#Object\ttarget\tdist" >  $TMP.tsv
cat $TMP.dist                   >> $TMP.tsv

section "tsv_group"
$THIS/../tsv/tsv_group $TMP.tsv  -by "Object"  -min "dist"  > $TMP.gr
  # Object dist_min

section "report"
$THIS/../tsv/tsv_rename $TMP.gr 2 "dist" -qc > $TMP.gr1
  # Object dist
$THIS/../tsv/tsv_expand.sh $TMP.gr1 $TMP.tsv ''
  # Object dist target

cut -f 1,2 "Metadata.tsv" > $TMP.meta
  # Object accver
$THIS/../tsv/tsv_expand.sh $TMP.gr1 $TMP.meta ''
  # Object dist target accver
cut -f 1 --complement $TMP.gr1 > $TMP.gr2
  # dist target accver
$THIS/../tsv/tsv_rename $TMP.meta  1 "target"       -qc > $TMP.meta1
$THIS/../tsv/tsv_rename $TMP.meta1 2 "$TARGET_NAME" -qc > $TMP.meta2
  # target "$TARGET_NAME"
$THIS/../tsv/tsv_expand.sh $TMP.gr2 $TMP.meta2 ''
  # dist target accver "$TARGET_NAME"
awk -F '\t' '{OFS="\t"; print $3, $4, $1};' $TMP.gr2 > $OUT


rm $TMP*
