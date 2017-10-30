#!/bin/csh -f

if ($# != 2) then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  echo "#2: tree version"
  exit 1
endif



set Leaves = `grep '^# Leaves: ' $1/old/makeDistTree.$2`
if ($?) exit 1

set CHRON_Leaves = `grep '^CHRON: Optimizing new leaves: ' $1/old/makeDistTree.$2`
if ($?) exit 1

set CHRON_topology = `grep '^CHRON: Initial topology: ' $1/old/makeDistTree.$2`
if ($?) exit 1

set Radius = `grep '^Ave. radius: ' $1/old/makeDistTree.$2`
if ($?) exit 1

echo $Leaves[3] $CHRON_Leaves[5] $CHRON_topology[4] $Radius[2]
