#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Compare the criteria; exit 1 if they are different"
  echo "#1: makeDistTree output #1 ('NEW')"
  echo "#2: makeDistTree output #2 ('OLD')"
  exit 1
fi
OUT1=$1
OUT2=$2


L1=(`grep "^OUTPUT" -A 1  $OUT1 | tail -1`)
L2=(`grep "^OUTPUT" -A 1  $OUT2 | tail -1`)
A=${L1[7]}
B=${L2[7]}
if [ $A == "0.000" -a $B == "0.000" ]; then
  exit 0
fi
R=`echo "($A - $B) / $A * 1000" | bc -l | cut -c1`
if [ $R == '-' ]; then
  R=`echo "($B - $A) / $A * 1000" | bc -l | cut -c1`
fi
# | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g'
if [ $R != '.' -a $R != '0' ]; then
  echo "NEW: ${L1[*]}"
  echo "OLD: ${L2[*]}"
  exit 1
fi

