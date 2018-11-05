#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Process #1/{hybrid,hybrid_parent_pairs}"
  echo "#1: incremental distance tree directory"
  exit 1
fi

INC=$1


VER=`cat $INC/version`


if [ -e $INC/hybrid_parent_pairs ]; then
  mv $INC/hybrid_parent_pairs $INC/hist/hybrid_parent_pairs.$VER
fi


if [ ! -e $INC/hybrid ]; then
  exit 0
fi
wc -l $INC/hybrid

if [ ! -s $INC/hybrid ]; then
  rm $INC/hybrid
  exit 0
fi

$INC/hybrid2db.sh $INC/hybrid

cat $INC/hybrid | awk '$7 == 1' | cut -f 1 >  $INC/outlier.add
cat $INC/hybrid | awk '$8 == 1' | cut -f 3 >> $INC/outlier.add
cat $INC/hybrid | awk '$9 == 1' | cut -f 4 >> $INC/outlier.add
uniq.sh $INC/outlier.add

$INC/objects_in_tree.sh $INC/outlier.add 0
trav -noprogress $INC/outlier.add "cp /dev/null $INC/outlier/%f"
rm $INC/outlier.add

mv $INC/hybrid $INC/hist/hybrid.$VER
