#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print a random sample of a .dm-file "
  echo "#1: Input .dm file without .dm"
  echo "#2: Sample seed"
  echo "#3: Sample size"
  exit 1
fi
INPUT=$1
SEED=$2
SIZE=$3


TMP=`mktemp`


$THIS/dm2objs $INPUT > $TMP.list
$THIS/../setRandOrd $TMP.list  -seed $SEED  -sigpipe | head -$SIZE | sort > $TMP.subset
$THIS/dm2subset $INPUT $TMP.subset


rm $TMP*
