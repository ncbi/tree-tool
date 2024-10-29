#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 3 ]; then
  echo "#1: List1"
  echo "#2: List2"
  echo "#3: number (0/1)"
  echo "Print set-theoretic intersection of List1 and List2"
  exit 1
fi
L1=$1
L2=$2
NUMP=$3


NUM=""
if [ $NUMP == 1 ]; then
  NUM="-number"
fi

TMP=$( mktemp )
#comment $TMP
#set -x

$THIS/setMinus $NUM $L1 $L2 > $TMP
$THIS/setMinus $NUM $L1 $TMP

rm $TMP*




