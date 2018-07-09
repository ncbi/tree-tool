#!/bin/csh -f

if ($# != 2) then
  echo "Initialize an incremental distance tree"
  echo "#1: Output directory"
  echo "#2: grid_min (> 0)"
  exit 1
endif



mkdir $1
if ($?) exit 1


cp /dev/null $1/tree
if ($?) exit 1

cp /dev/null $1/dissim
if ($?) exit 1

cp /dev/null $1/leaf
if ($?) exit 1

mkdir $1/new
if ($?) exit 1

cp /dev/null $1/delete-hybrid
if ($?) exit 1

mkdir $1/search
if ($?) exit 1

mkdir $1/outlier
if ($?) exit 1

echo "1" > $1/version
if ($?) exit 1

mkdir $1/hist
if ($?) exit 1

echo $2 > $1/grid_min
if ($?) exit 1

cp /dev/null $1/runlog
if ($?) exit 1

echo "exit 1" > $1/request2dissim.sh
if ($?) exit 1
chmod a+x $1/request2dissim.sh
if ($?) exit 1

echo "exit 1" > $1/objects_in_tree.sh
if ($?) exit 1
chmod a+x $1/objects_in_tree.sh
if ($?) exit 1

echo "exit 1" > $1/request_closest.sh
if ($?) exit 1
chmod a+x $1/request_closest.sh
if ($?) exit 1

echo "exit 1" > $1/hybrid2db.sh
if ($?) exit 1
chmod a+x $1/hybrid2db.sh
if ($?) exit 1

echo "exit 1" > $1/db2hybrid.sh
if ($?) exit 1
chmod a+x $1/db2hybrid.sh
if ($?) exit 1

echo "exit 1" > $1/db2unhybrid.sh
if ($?) exit 1
chmod a+x $1/db2unhybrid.sh
if ($?) exit 1
