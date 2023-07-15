#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Create a phylogenetic tree from a FASTA files"
  echo "#1: FASTA"
  echo "#2: #1 is proteins (0/1)"
  echo "#3: DNA strand is known (0/1)"
  echo "#4: global alignment (0/1)"
  echo "#5: `fasta2dissim -help | grep ^-distance -A 1 | tail -1 | sed 's/^ *//1'`"
  echo "#6: output tree"
  exit 1
fi
FASTA=$1
PROT=$2
KNOWN_STRAND=$3
GLOB=$4
DIST=$5
TREE=$6


TMP=`mktemp`
comment $TMP 
#set -x


section "Computing dissimilarities"

PAR="-distance $DIST"
if [ $GLOB == 1 ]; then
  PAR="$PAR  -global"
fi
if [ $PROT == 0 -a $KNOWN_STRAND == 0 ]; then
  PAR="$PAR  -unknown_strand"
fi

if [ $PROT == 1 ]; then
  POWER=""
  if [ $DIST != "diff" ]; then
    POWER="-power 0.5"
  fi
  $THIS/../dissim/fasta2dissim  -threads 30  $FASTA  $PAR  -aa  -blosum62  $POWER  -dataset $TMP  
else
  $THIS/../dissim/fasta2dissim  -threads 30  $FASTA  $PAR  -dataset $TMP  
fi


section "Building tree"
# PAR
$THIS/makeDistTree  -data $TMP  -dissim_attr "dissim"  -variance pow  -variance_power 5  -optimize  -subgraph_iter_max 5  -output_tree $TREE  -threads 10


rm $TMP*

