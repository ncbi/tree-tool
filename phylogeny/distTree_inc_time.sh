#!/bin/csh -f

if ($# != 2) then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  echo "#2: tree version"
  exit 1
endif



set L1 = `grep '^# Leaves: ' $1/old/makeDistTree.$2`
if ($?) exit 1

set L2 = `grep '^CHRON: Optimizing new leaves: ' $1/old/makeDistTree.$2`
if ($?) exit 1

set L3 = `grep '^CHRON: Initial topology: ' $1/old/makeDistTree.$2`
if ($?) exit 1

echo $L1[3] $L2[5] $L3[4]
