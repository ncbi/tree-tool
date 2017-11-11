#!/bin/csh -f

if ($# != 3) then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  echo "#2: tree version"
  echo "#3: Chronometer::enabled"
  exit 1
endif



set Leaves = `grep '^# Leaves: ' $1/old/makeDistTree.$2`
if ($?) exit 1

set Dissim = `grep '^# Dissimilarities: ' $1/old/makeDistTree.$2`
if ($?) exit 1

set CHRON_topology = `grep '^CHRON: Initial topology: ' $1/old/makeDistTree.$2`
if ($?) exit 1

set CHRON_Leaves = `grep '^CHRON: Optimizing new leaves: ' $1/old/makeDistTree.$2`
if ($?) exit 1


set CHRON_tree2subgraph = 0
set CHRON_tree2subgraphDissim = 0
set CHRON_subgraph = 0
set CHRON_subgraph2tree = 0
if ($3) then
  set L = `grep '^CHRON: Tree -> subgraph: ' $1/old/makeDistTree.$2`
  if ($?) exit 1
  set CHRON_tree2subgraph = $L[5]
  
  set L = `grep '^CHRON: Tree -> subgraph dissim: ' $1/old/makeDistTree.$2`
  if ($?) exit 1
  set CHRON_tree2subgraphDissim = $L[6]
  
  set L = `grep '^CHRON: Subgraph optimization: ' $1/old/makeDistTree.$2`
  if ($?) exit 1
  set CHRON_subgraph = $L[4]
  
  set L = `grep '^CHRON: Subgraph -> tree: ' $1/old/makeDistTree.$2`
  if ($?) exit 1
  set CHRON_subgraph2tree = $L[5]
endif


set Radius = `grep '^Ave. radius: ' $1/old/makeDistTree.$2`
if ($?) exit 1


#     #Leaves   #Dissim           topology           Leaves           tree2subgraph        tree2subgraphDissim        subgraph        subgraph2tree  Radius   
echo $Leaves[3] $Dissim[3] $CHRON_topology[4] $CHRON_Leaves[5] $CHRON_tree2subgraph $CHRON_tree2subgraphDissim $CHRON_subgraph $CHRON_subgraph2tree $Radius[3]
