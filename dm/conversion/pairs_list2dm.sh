#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print a .dm-file out of a list of two-way attribute files"
  echo "File line format: <obj1> <obj2> <value> [<comment>]"
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


$THIS/../trav $LIST "cat %f | tr ' ' '\t' | cut -f 1" >  $TMP.obj_raw
$THIS/../trav $LIST "cat %f | tr ' ' '\t' | cut -f 2" >> $TMP.obj_raw
sort -u $TMP.obj_raw > $TMP.obj

N=`cat $TMP.obj | wc -l`
if [ $N == 0 ]; then
  error "No data"
fi

echo "ObjNum $N name nomult"
echo "Attributes"

ATTR="Real"
if [ $POS == 1 ]; then
  ATTR="Positive2"
fi
ATTR="$ATTR $DEC"
$THIS/../trav $LIST "basename %f" | sed 's/$/ '"$ATTR"'/1' | sed 's/^/  /1'

echo "DATA"
cat $TMP.obj

$THIS/../trav $LIST "basename %f; echo PAIRS; wc -l %f | awk %q{print %D1};%q; sort %f"


rm $TMP*  
