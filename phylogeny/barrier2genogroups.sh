#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
source $THIS/../qsub_env.sh
COEFF=100
if [ $# -ne 5 ]; then
  echo "Print dependence of genogroup size on genogroup barrier"
  echo "#1: distance tree"
  echo "#2: min. barrier * $COEFF (integer)"
  echo "#3: max. barrier * $COEFF (integer)" 
  echo "#4: step"
  echo "#5: output directory"
  exit 1
fi
TREE=$1
MIN=$2
MAX=$3
STEP=$4
OUT=$5


$THIS/../grid_wait.sh 1 >> /dev/stderr
i=$MIN
while [ $i -le $MAX ]; do
  T=$( echo "$i/$COEFF" | bc -l | sed 's/0*$//1' | sed 's/^\./0./1' )
  $QSUB_5 -N j$i "$THIS/tree2genogroup $TREE $T  -genogroup_table $OUT/$T" > /dev/null
  i=$(( i + STEP ))
done
$THIS/../qstat_wait.sh 2000 1

$THIS/../trav -noprogress $OUT 'echo "%f %gcut -f 2 %d/%f | sort -u | wc -l%g"'
