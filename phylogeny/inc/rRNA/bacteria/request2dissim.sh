#!/bin/bash
source CPP_DIR/bash_common.sh
if [ $# -ne 3 ]; then
  echo "$0"
  echo "#1: dissimilaruty requests (input)"
  echo "#2: dissim (output)"
  echo "#3: log"
  exit 1
fi
REQUEST=$1
DISSIM=$2
LOG=$3


INC=`dirname $0`
# PAR
CPP_DIR/dissim/dna_pair2dissim  -log $LOG  -coeff 0.00155  $REQUEST $INC/../seq 600 $DISSIM
  # was: 1200
rm -f $LOG
