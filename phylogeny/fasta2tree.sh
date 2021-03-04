#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Create a phylogenetic tree from a FASTA files"
  echo "#1: FASTA"
  echo "#2: #1 is proteins (1/0)"
  echo "#3: output tree"
  exit 1
fi
FASTA=$1
PROT=$2
TREE=$3


TMP=`mktemp`


section "Computing dissimilarities ..."
# PAR
if [ $PROT == 1 ]; then
  $THIS/../dissim/fasta2dissim  $FASTA  -aa  -global  -blosum62  -power 0.5  -dataset $TMP  
else
  error "Not implemented for DNA"
fi

echo ""
section "Builing tree ..."
# PAR
$THIS/makeDistTree  -data $TMP  -dissim_attr "dissim"  -variance "linExp"  -optimize  -subgraph_iter_max 5  -output_tree $TREE


rm $TMP*

