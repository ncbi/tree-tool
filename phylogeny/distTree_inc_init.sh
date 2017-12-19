#!/bin/csh -f

if ($# != 1) then
  echo "Initialize an incremental distance tree"
  echo "#1: Output distance tree data directory"
  exit 1
endif


mkdir $1
if ($?) exit 1

cp /dev/null $1/dissim
if ($?) exit 1

cp /dev/null $1/leaf
if ($?) exit 1

cp /dev/null $1/tree
if ($?) exit 1

mkdir $1/outlier
if ($?) exit 1

#cp /dev/null $1/strong_outliers
#if ($?) exit 1

mkdir $1/new
if ($?) exit 1

mkdir $1/search
if ($?) exit 1

mkdir $1/old
if ($?) exit 1

echo "1" >> $1/version
if ($?) exit 1

# To be modified
echo "exit 1" > $1/request2dissim.sh
if ($?) exit 1
chmod a+x $1/request2dissim.sh
if ($?) exit 1
