#!/bin/csh -f

if ($# != 3) then
  echo "Process #1/hybrid"
  echo "Output: #1/delete-hybrid"
  echo "#1: incremental distance tree directory"
  echo "#2: tree version"
  echo "#3: List of objects to test hybridness against #1/hybrid, sorted"
  exit 1
endif


if (! -z $1/delete-hybrid)  exit 1

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
mv $1/hybrid $1/hist/hybrid.$2
if ($?) exit 1

setMinus $3 $1/outlier.add > $1/test.list
if ($?) exit 1

rm $1/outlier.add

# hybrid-add
set HYBRIDNESS_MIN = `cat $1/hybridness_min`
mkdir $1/hybrid
if ($?) exit 1
mkdir $1/log
if ($?) exit 1
# grid ??
trav $1/test.list "distTree_inc_new2hybrid.sh %f $1 $HYBRIDNESS_MIN $1/hybrid/%f $1/log/%f"
if ($?) exit 1
rmdir $1/log
if ($?) exit 1
trav -noprogress $1/hybrid "cat $1/hybrid/%f" > $1/hybrid-add
if ($?) exit 1
rm -r $1/hybrid/
if ($?) exit 1

rm $1/test.list

# Cf. processing $1/hybrid
wc -l $1/hybrid-add
$1/hybrid2db.sh $1/hybrid-add
if ($?) exit 1
cut -f 1 $1/hybrid-add | sort > $1/delete-hybrid
if ($?) exit 1
$1/objects_in_tree.sh $1/delete-hybrid 0
if ($?) exit 1
trav -noprogress $1/delete-hybrid "cp /dev/null $1/outlier/%f"
if ($?) exit 1
mv $1/hybrid-add $1/hist/hybrid-add.$2
if ($?) exit 1







