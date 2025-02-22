#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Create a phylogenetic tree from a FASTA files"
  echo "Output: #1.tree, #1.dm"
  echo "#1: FASTA where identifiers start with CLASS-"
  echo "#2: output .qual-file"
  exit 1
fi
F=$1
OUT=$2


TMP=$( mktemp )
comment $TMP


$THIS/fasta2tree.sh $F 1 0 1 "min_edit" $TMP
  # should have been semiglobal ??
  
super_section "Quality"
mkdir $TMP.phen
$THIS/../genetics/fa2list.sh $F > $TMP.list
trav $TMP.list "echo %f | cut -f 1 -d '-' | sed 's/$/ 0/1' > $TMP.phen/%f" -threads 10
$THIS/tree_quality_phen.sh $TMP.tree  '' $TMP.phen 0 1 $OUT

echo "Bad ortholog classes (may have outliers):"
awk '$4 != 0' $OUT


rm -r $TMP*
