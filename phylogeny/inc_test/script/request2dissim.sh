#!/bin/bash
THIS=`dirname $0`
source bash_common.sh
if [ $# -ne 3 ]; then
  echo "Compute dissimilarities"
  echo "#1: file with pairs <Object1> <Object2>"
  echo "#2: output file"
  echo "#3: log (temporary)"
  exit 1
fi
REQ=$1
OUT=$2
LOG=$3


DT/phylogeny/hash_request2dissim $REQ $THIS/../hash $OUT  -intersection_min 50  -ratio_min 0.5  -log $LOG 
