#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Create a phylogenetic tree from a FASTA files"
  echo "Output: #1.tree, #1.dm"
  echo "#1: FASTA"
  echo "#2: #1 is proteins (0/1)"
  echo "#3: DNA strand is known (0/1)"
  echo "#4: global alignment (0/1)"
  echo "#5: `$THIS/../dissim/fasta2dissim -help | grep ^-distance -A 1 | tail -1 | sed 's/^ *//1'`"
  echo "#6: output file prefix"
  exit 1
fi
FASTA=$1
PROT=$2
KNOWN_STRAND=$3
GLOB=$4
DIST=$5
OUT=$6


section "Computing dissimilarities"

PAR="-distance $DIST"
if [ $GLOB == 1 ]; then
  PAR="$PAR  -global"
fi
if [ $PROT == 1 ]; then
  POWER=""  # ??
  if [ $DIST != "diff" ]; then
    POWER="-power 0.5"
  fi
  PAR="$PAR  -aa  -blosum62"
else
  if [ $KNOWN_STRAND == 0 ]; then
    PAR="$PAR  -unknown_strand"
  fi
fi

$THIS/../dissim/fasta2dissim  -threads 30  $FASTA  $PAR  -dataset $OUT.dm


section "Building tree"
# PAR
$THIS/makeDistTree  -data $OUT.dm  -dissim_attr "dissim"  -variance pow  -variance_power 5  -optimize  -subgraph_iter_max 5  -output_tree $OUT.tree  -threads 10


