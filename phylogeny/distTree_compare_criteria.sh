#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Compare the criteria; exit 1 if they are different"
  echo "#1: makeDistTree output #1"
  echo "#2: makeDistTree output #2"
  exit 1
fi
OUT1=$1
OUT2=$2


L1=`grep "^OUTPUT" -A 1  $OUT1 | tail -1`
L2=`grep "^OUTPUT" -A 1  $OUT2 | tail -1`
if [ "$L1" -ne "$L2" ]; then
  echo "$L1"
  echo "$L2"
  exit 1
fi

