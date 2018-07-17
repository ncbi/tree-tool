#!/bin/csh -f

if ($# != 2) then
  echo "Process #1/hybrid"
  echo "#1: incremental distance tree directory"
  echo "#2: tree version"
  exit 1
endif


if (-e $1/unhybrid)  exit 1


$1/db2unhybrid.sh $1/unhybrid
if ($?) exit 1

if (-z $1/unhybrid) then
  rm $1/unhybrid
else
  echo ""
  wc -l $1/unhybrid
  trav $1/unhybrid "mv $1/outlier/%f $1/new/"
	if ($?) exit 1
  $1/objects_in_tree.sh $1/unhybrid null
	if ($?) exit 1
	mv $1/unhybrid $1/hist/unhybrid.$2
	if ($?) exit 1
endif

