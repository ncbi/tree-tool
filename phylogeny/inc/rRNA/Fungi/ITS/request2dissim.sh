#!/bin/bash
source CPP_DIR/bash_common.sh
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


INC=`dirname $0`
CPP_DIR/dissim/dna_pair2dissim  -log $LOG  -coeff 0.0082  $REQUEST $INC/../seq 140 $DISSIM
#dna_pairs2dissim  -log $LOG  -power 1.01  $REQUEST /home/brovervv/panfs/marker/Fungi/ITS/seq 140 $DISSIM
  # Needs recomputation of inc/dissim
rm -f $LOG

