#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
META="Metadata.tsv"
if [ $# -ne 2 ]; then
  echo "Print match of a genogroup barrier with taxonomic species"
  echo "Input: $META"
  echo "#1: distance tree"
  echo "#2: species barrier"
  exit 1
fi
TREE=$1
BARR=$2


$THIS/../check_file.sh $META 1


TMP=$( mktemp )
#comment $TMP


cut -f 1,7 $META | awk '$2 != ""' > $TMP.species

echo -e "#Object\tgenogroup" > $TMP.gr
$THIS/tree2genogroup $TREE $BARR  -genogroup_table $TMP
cat $TMP >> $TMP.gr

$THIS/../tsv/tsv_expand.sh $TMP.gr $TMP.species '' &> /dev/null
cut -f 2,3 $TMP.gr | sort | uniq -c | awk -F ' ' '{OFS="\t"; print $2, $3, $1};' > $TMP.uniq
$THIS/../dm/assignment $TMP.uniq -qc


rm $TMP* 

