#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Hash difference"
  echo "#1: file with hashes"
  echo "#2: file with hashes"
  exit 1
fi
F1=$1
F2=$2


sort -ncu $F1
sort -ncu $F2


N1=`cat $F1 | wc -l`
N2=`cat $F2 | wc -l`
COMMON=`$THIS/../setIntersect.sh $F1 $F2 1 | wc -l`
DIFF=$(( $N1 - $COMMON + $N2 - $COMMON ))
echo -e "$F1\t$F2\t$N1\t$N2\t$DIFF"

