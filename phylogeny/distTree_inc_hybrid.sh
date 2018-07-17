#!/bin/csh -f

if ($# != 2) then
  echo "Process #1/hybrid"
  echo "#1: incremental distance tree directory"
  echo "#2: tree version"
  exit 1
endif


if (! -e $1/hybrid) exit 0

if (-z $1/hybrid) then
  rm $1/hybrid
  exit 0
endif

wc -l $1/hybrid
$1/hybrid2db.sh $1/hybrid
if ($?) exit 1
cut -f 1 $1/hybrid | sort > $1/outlier.add
if ($?) exit 1
$1/objects_in_tree.sh $1/outlier.add 0
if ($?) exit 1
trav -noprogress $1/outlier.add "cp /dev/null $1/outlier/%f"
if ($?) exit 1
rm $1/outlier.add
if ($?) exit 1
mv $1/hybrid $1/hist/hybrid.$2
if ($?) exit 1
