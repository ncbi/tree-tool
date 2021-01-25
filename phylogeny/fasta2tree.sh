#!/bin/bash
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


# PAR
if [ $PROT == 1 ]; then
  $THIS/../dissim/fasta2dissim  -aa  -global  -blosum62  -power 0.5  -dataset $TMP  $FASTA
else
  error "Not implemented yet for DNA"
fi

$THIS/makeDistTree  -data $TMP  -dissim_attr "dissim"  -optimize  -output_tree $TREE


rm -f $TMP*

