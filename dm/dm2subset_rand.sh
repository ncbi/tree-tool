#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Seed sample of a .dm-file "
  echo "#1: Input .dm file without .dm"
  echo "#2: Sample seed"
  echo "#3: Sample size"
  exit 1
fi
INPUT=$1
BASE_SEED=$2
BASE_SIZE=$3


BASE=$INPUT-$BASE_SEED
$THIS/dm2objs $INPUT | sort > $INPUT.list
$THIS/../setRandOrd $INPUT.list  -seed $BASE_SEED | head -$BASE_SIZE | sort > $BASE.list
$THIS/dm2subset $INPUT $BASE.list > $BASE.dm
rm $BASE.list


