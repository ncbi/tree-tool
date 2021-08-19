#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Create a phylogenetic tree from a FASTA files"
  echo "#1: FASTA"
  echo "#2: #1 is proteins (1/0)"
  echo "#3: global alignment (0/1)"
  echo "#4: output tree"
  exit 1
fi
FASTA=$1
PROT=$2
GLOB=$3
TREE=$4


TMP=`mktemp`


section "Computing dissimilarities"
GLOB_PAR=""
if [ $GLOB == 1 ]; then
  GLOB_PAR="-global"
fi
# PAR
if [ $PROT == 1 ]; then
  $THIS/../dissim/fasta2dissim  $FASTA  $GLOB_PAR  -aa  -blosum62  -power 0.5  -dataset $TMP  
else
  $THIS/../dissim/fasta2dissim  $FASTA  $GLOB_PAR  -dataset $TMP  
fi

section "Builing tree"
# PAR
$THIS/makeDistTree  -data $TMP  -dissim_attr "dissim"  -variance "linExp"  -optimize  -subgraph_iter_max 5  -output_tree $TREE


rm $TMP*

