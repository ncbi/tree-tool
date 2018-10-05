#!/bin/bash
set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C

if [ $# != 2 ]; then
  echo "Process #1/{hybrid,hybrid_parent_pairs}"
  echo "#1: incremental distance tree directory"
  echo "#2: tree version"
  exit 1
fi

inc=$1
ver=$2


if [ -e $inc/hybrid_parent_pairs ]; then
  mv $inc/hybrid_parent_pairs $inc/hist/hybrid_parent_pairs.$ver
fi


if [ ! -e $inc/hybrid ]; then
  exit 0
fi
wc -l $inc/hybrid

if [ ! -s $inc/hybrid ]; then
  rm $inc/hybrid
  exit 0
fi

$inc/hybrid2db.sh $inc/hybrid

cat $inc/hybrid | awk '$7 == 1' | cut -f 1 >  $inc/outlier.add
cat $inc/hybrid | awk '$8 == 1' | cut -f 3 >> $inc/outlier.add
cat $inc/hybrid | awk '$9 == 1' | cut -f 4 >> $inc/outlier.add
uniq.sh $inc/outlier.add

$inc/objects_in_tree.sh $inc/outlier.add 0
trav -noprogress $inc/outlier.add "cp /dev/null $inc/outlier/%f"
rm $inc/outlier.add

mv $inc/hybrid $inc/hist/hybrid.$ver
