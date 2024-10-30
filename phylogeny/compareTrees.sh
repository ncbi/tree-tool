#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Robinson-Foulds comparison"
  echo "#1: tree 1"
  echo "#2: tree 2"
  exit 1
fi
T1=$1
T2=$2


TMP=$( mktemp )
#comment $TMP
#set -x


$THIS/compareTrees $T1 $T2 -qc > $TMP.out
grep -w '^match+' $TMP.out > $TMP.1 || true
grep -w '^match-' $TMP.out > $TMP.2 || true

M0=$( < $TMP.1 wc -l )
M1=$( < $TMP.2 wc -l )
#echo $M0
#echo $M1

P=$( echo "scale=2; $M1 * 100 / ($M0 + $M1)" | bc -l )
echo -e "$T1\t$T2\t$M1\t$M0\t$P"


rm $TMP*
