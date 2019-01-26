#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "#1: Tree 1"
  echo "#2: Tree 2"
  echo "#3: Collapse rare branches: directed|undirected|none"
  exit 1
fi


TMP=`mktemp`


$THIS/compareTrees $1 $2  -frequency $3 > $TMP.out

M=(`grep -w match $TMP.out | tr '\t' ' '  | sed 's/ .*$//1' | sort | uniq -c`)
#   5177 match+
#   2092 match-
M0=${M[0]}
M1=${M[2]}

N=$(( $M0 + $M1 ))
P=`echo "scale=2; $M1 * 100 / ($M0 + $M1)" | bc -l`
echo "$1	$2	$M1 / ($M0 + $M1) = $P %"


rm -f $TMP*
