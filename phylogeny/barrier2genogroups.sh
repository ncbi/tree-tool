#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Print dependence of genogroup size on genogroup barrier"
  echo "#1: distance tree"
  echo "#2: min. barrier"
  echo "#3: max. barrier" 
  echo "#4: step"
  echo "#5: output directory"
  echo "#6: output statistics: dist_max num_clusters"
  exit 1
fi
TREE=$1
MIN=$2
MAX=$3
STEP=$4
DIR=$5
OUT=$6


T=$MIN
while true
do
  R=$( echo "$T > $MAX" | bc )
  if [ $R == 1 ]; then
    break
  fi
  printf "\r$T" > /dev/stderr
  $THIS/tree2genogroup $TREE $T  -genogroup_table $DIR/$T > /dev/null 
  T=$( echo "$T + $STEP" | bc -l | sed 's/^\./0./1' )
done
echo "" > /dev/stderr

$THIS/../trav $DIR 'echo -e "%f\t%D( cut -f 2 %d/%f | sort -u | wc -l )"' > $OUT
