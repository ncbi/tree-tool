#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Robinson-Foulds comparison"
  echo "#1: Tree 1"
  echo "#2: Tree 2"
  exit 1
fi
T1=$1
T2=$2


TMP=`mktemp`
#echo $TMP
#set -x


$THIS/compareTrees $T1 $T2 > $TMP.out

M=(`grep -w 'match' $TMP.out | tr '\t' ' '  | sed 's/ .*$//1' | sort | uniq -c`)
#   5177 match+
#   2092 match-
M0=${M[0]}
M1=${M[2]}

#head -20 $TMP.out

P=`echo "scale=2; $M1 * 100 / ($M0 + $M1)" | bc -l`
echo "$1	$2	$M1 / ($M0 + $M1) = $P %"


rm $TMP*
