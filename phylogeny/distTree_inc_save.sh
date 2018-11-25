#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Save the parameters and scripts of an incremental distance tree"
  echo "#1: Incremental tree directory"
  echo "#2: Output directory"
  exit 1
fi
INC=$1
OUT=$2


rm -rf $OUT
mkdir $OUT
for F in `ls $INC/`; do
  if [ -d $INC/$F ]; then
    continue;
  fi

  set +o errexit
  N=`file $INC/$F | grep -c ': data$'`
  set -o errexit
  if [ $N -ne 0 ]; then
    continue
  fi
  
  N=`cat $INC/$F | wc -l`
  if [ $N -lt 500 ]; then  # PAR
    cp $INC/$F $OUT/
  fi
done
