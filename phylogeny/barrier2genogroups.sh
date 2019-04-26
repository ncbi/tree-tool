#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Print dependence of genogroup size on genogroup barrier"
  echo "#1: incremental distance tree directory"
  echo "#2: min. barrier * 10 (integer)"
  echo "#3: max. barrier * 10 (integer)" 
  echo "#4: output directory"
  exit 1
fi
INC=$1
MIN=$2
MAX=$3
OUT=$4


i=$MIN
while [ $i -le $MAX ]; do
  T=`echo "$i/10" | bc -l | sed 's/0*$//1' | sed 's/^\./0./1'`
  $QSUB_5 -N j$i "$THIS/tree2genogroup $INC/tree $T  -genogroup_table $OUT/$T" > /dev/null
  i=$(( $i + 1 ))
done
$THIS/../qstat_wait.sh 2000 1

$THIS/../trav -noprogress $OUT 'echo "%f `cut -f 2 %d/%f | sort | uniq | wc -l`"'
