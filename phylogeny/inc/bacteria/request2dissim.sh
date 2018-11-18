#!/bin/bash
source bash_common.sh
if [ $# -ne 3 ]; then
  echo "Compute dissimilarities"
  echo "#1: file with pairs <Object1> <Object2>"
  echo "#2: Output file"
  echo "#3: log (temporary)"
  exit 1
fi

DIR=/home/brovervv/panfs/GenBank/bacteria
# PAR
~/code/genetics/combine_dissims3.sh $1 $2 $DIR/hash 50 0.5 $DIR/dissim_scale $DIR/hmm-univ.stat 0 0.57 $3


