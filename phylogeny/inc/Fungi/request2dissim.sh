#!/bin/csh -f

if ($# != 3) then
  echo "Compute dissimilarities"
  echo "#1: file with pairs <Object1> <Object2>"
  echo "#2: Output file"
  echo "#3: log (temporary)"
  exit 1
endif

set DIR = /home/brovervv/panfs/GenBank/Fungi
# was:                                             200
~/code/genetics/combine_dissims.sh $1 $2 $DIR/hash  10 0.1 $DIR/dissim_scale $DIR/hmm-univ.stat 1 0.57 $3
#                                  1  2  3          4  5   6                 7                  8 9    10
if ($?) exit 1



