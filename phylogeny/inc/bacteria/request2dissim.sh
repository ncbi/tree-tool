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


DIR=/home/brovervv/panfs/GenBank/bacteria
# PAR
~/code/genetics/combine_dissims3.sh $REQ $OUT $DIR/hash 50 0.5 $THIS/dissim_scale $THIS/hmm-univ.stat 0 0.57 $LOG


