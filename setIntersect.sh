#!/bin/bash
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 3 ]; then
  echo "#1 - List1"
  echo "#2 - List2"
  echo "#3 - number(0/1)"
  echo "Make set-theoretic intersection of List1 and List2"
  exit 1
fi


NUM=""
if [ $3 == 1 ]; then
  NUM="-number"
fi

TMP=`mktemp` 
#echo $TMP

$THIS/setMinus $NUM $1 $2 > $TMP
$THIS/setMinus $NUM $1 $TMP

rm $TMP*




