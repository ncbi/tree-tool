#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print a .dm-file out of a list square matrices each containing a subset of objects"
  echo "Format of a matrix: <number of objects in subset> \n {<object name> <number>+ \n}+"
  echo "File names become attribute names"
  echo "#1: List of matrix files"
  echo "#2: Numbers are non-negative (0/1)"
  echo "#3: Decimals"
  exit 1
fi
LIST=$1
POS=$2
DEC=$3


TMP=`mktemp`
#echo $TMP  


trav $LIST "tail -n +2 %f" | tr ' ' '\t' | cut -f1 | sort | uniq > $TMP.obj

N=`cat $TMP.obj | wc -l`
if [ $N == 0 ]; then
  echo "No data"
  exit 1
fi

echo "ObjNum $N name nomult"
echo "Attributes"

ATTR="Real"
if [ $POS == 1 ]; then
  ATTR="Positive2"
fi
ATTR="$ATTR $DEC"
trav $LIST "basename %f" | sed 's/$/ '"$ATTR"'/1' | sed 's/^/  /1'

echo "DATA"
cat $TMP.obj

trav $LIST "basename %f; echo PARTIAL; cat %f"


rm $TMP*  
