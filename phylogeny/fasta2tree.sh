#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Create a phylogenetic tree from a FASTA files"
  echo "#1: FASTA"
  echo "#2: #1 is proteins (1/0)"
  echo "#3: DNA strand is known (0/1)"
  echo "#4: global alignment (0/1)"
  echo "#5: output tree"
  exit 1
fi
FASTA=$1
PROT=$2
KNOWN_STRAND=$3
GLOB=$4
TREE=$5


TMP=`mktemp`


section "Computing dissimilarities"

PAR=""
if [ $GLOB == 1 ]; then
  PAR="-global"
fi
if [ $PROT == 0 -a $KNOWN_STRAND == 0 ]; then
  PAR="$PAR -unknown_strand"
fi

# PAR
VAR="linExp"
if [ $PROT == 1 ]; then
  $THIS/../dissim/fasta2dissim  -threads 30  $FASTA  $PAR  -aa  -blosum62  -power 0.5  -dataset $TMP  
else
  $THIS/../dissim/fasta2dissim  -threads 30  $FASTA  $PAR  -dataset $TMP  
  VAR="pow  -variance_power 2"
fi


section "Building tree"
# PAR
$THIS/makeDistTree  -data $TMP  -dissim_attr "dissim"  -variance $VAR  -optimize  -subgraph_iter_max 5  -output_tree $TREE  -threads 10


rm $TMP*

