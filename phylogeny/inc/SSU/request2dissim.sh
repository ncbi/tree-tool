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

dna_pairs2dissim  -log $LOG  -coeff 25  -relative  $REQUEST /home/brovervv/panfs/marker/SSU/seq 300 $DISSIM
  # 0.01 * 2500 = 25
rm -f $LOG

