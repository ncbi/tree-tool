#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Create a phylogenetic tree from a FASTA file"
  echo "Output: #6.dm, #6.tree, #6.nw; for #2 = 1 and #5 = 'min_edit': #6.gr, #6.roots"
  echo "#1: FASTA"
  echo "#2: #1 is proteins (0/1)"
  echo "#3: DNA strand is known (0/1)"
  echo "#4: global alignment (0/1)"
  echo "#5: $( $THIS/../dissim/fasta2dissim -help | grep ^-distance -A 1 | tail -1 | sed 's/^ *//1' )"
  echo "#6: output file prefix"
  exit 1
fi
FASTA=$1
PROT=$2
KNOWN_STRAND=$3
GLOB=$4
DIST=$5
OUT=$6


super_section "Computing dissimilarities"

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

$THIS/../dissim/fasta2dissim  -threads 30  $FASTA  $PAR  -dataset $OUT


super_section "Building tree"
# PAR
$THIS/makeDistTree  -data $OUT  -dissim_attr "dissim"  -variance pow  -variance_power 5  -variance_dissim  -optimize  -subgraph_iter_max 5  -output_tree $OUT.tree  -threads 15
$THIS/printDistTree $OUT.tree  -order  -qc > $OUT.nw

if [ $DIST == "min_edit" ] && [ $PROT == 1 ] ; then
  LEN_AVG=$( $THIS/../genetics/fasta2len $FASTA -aa -qc | cut -f 2 | count | grep -w '^mean' | cut -f 2 )
  DIST_MAX=$( echo "$LEN_AVG / 370 * 4000" | bc -l )
    # PAR
    # For random protein sequences of length 370 aa the min_edit distance = ~4500
  echo "DIST_MAX=$DIST_MAX"
  $THIS/tree2genogroup $OUT.tree $DIST_MAX  -genogroup_table $OUT.gr
  cut -f 2 $OUT.gr | sort | uniq -c > $OUT.roots
fi


