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


set CHRON_getBestChange = 0
set CHRON_tree2subgraph = 0
set CHRON_tree2subgraphDissim = 0
set CHRON_subgraph = 0
set CHRON_subgraph2tree = 0
if ($3) then
  set L = `grep '^CHRON: getBestChange: ' $1/old/makeDistTree.$2`
  if ($?) exit 1
  set CHRON_getBestChange = $L[3]
  
  set L = `grep '^CHRON: tree2subgraph: ' $1/old/makeDistTree.$2`
  if ($?) exit 1
  set CHRON_tree2subgraph = $L[3]
  
  set L = `grep '^CHRON: tree2subgraphDissim: ' $1/old/makeDistTree.$2`
  if ($?) exit 1
  set CHRON_tree2subgraphDissim = $L[3]
  
  set L = `grep '^CHRON: subgraphOptimize: ' $1/old/makeDistTree.$2`
  if ($?) exit 1
  set CHRON_subgraph = $L[3]
  
  set L = `grep '^CHRON: subgraph2tree: ' $1/old/makeDistTree.$2`
  if ($?) exit 1
  set CHRON_subgraph2tree = $L[3]
endif


set Radius = `grep '^Ave. radius: ' $1/old/makeDistTree.$2`
if ($?) exit 1


#     #Leaves   #Dissim           topology           Leaves           getBestChange        tree2subgraph        tree2subgraphDissim        subgraph        subgraph2tree  Radius   
echo $Leaves[3] $Dissim[3] $CHRON_topology[4] $CHRON_Leaves[5] $CHRON_getBestChange $CHRON_tree2subgraph $CHRON_tree2subgraphDissim $CHRON_subgraph $CHRON_subgraph2tree $Radius[3]
