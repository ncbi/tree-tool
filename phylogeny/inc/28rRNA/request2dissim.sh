#!/bin/bash
source bash_common.sh
if [ $# -ne 3 ]; then
  echo "$0"
  echo "#1: dissimilarity requests (input)"
  echo "#2: dissim (output)"
  echo "#3: log"
  exit 1
fi
REQUEST=$1
DISSIM=$2
LOG=$3

dna_pairs2dissim  -log $LOG  -coeff 0.0066887  $REQUEST /home/brovervv/panfs/marker/Fungi/28SrRNA/seq 300 $DISSIM
  # was: 0.00634 
rm -f $LOG

