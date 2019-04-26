#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print time statistics for an incremental distance tree: # leaves, # dissimilarities, optimization time"
  echo "#1: Distance tree data"
  echo "#2: Tree version"
  echo "#3: Chronometer::enabled (0/1)"
  exit 1
fi
INC=$1
VER=$2
CHRON=$3


OUT=$INC/hist/makeDistTree.$VER
#echo $OUT


Leaves=(`grep '^# Leaves: ' $OUT`)
Dissim=(`grep '^# Dissimilarities: ' $OUT`)
CHRON_init=(`grep '^CHRON: Initial topology: ' $OUT`)
CHRON_optim=(`grep '^CHRON: Topology optimization: local: ' $OUT`)


CHRON_getBestChange=0
CHRON_tree2subgraph=0
CHRON_tree2subgraphDissim=0
CHRON_subgraph=0
CHRON_subgraph2tree=0
if [ $CHRON == 1 ]; then
  L=(`grep '^CHRON: getBestChange: ' $OUT`)
  CHRON_getBestChange=${L[2]}
  
  L=(`grep '^CHRON: tree2subgraph: ' $OUT`)
  CHRON_tree2subgraph=${L[2]}
  
  L=(`grep '^CHRON: tree2subgraphDissim: ' $OUT`)
  CHRON_tree2subgraphDissim=${L[2]}
  
  L=(`grep '^CHRON: subgraphOptimize: ' $OUT`)
  CHRON_subgraph=${L[2]}
  
  L=(`grep '^CHRON: subgraph2tree: ' $OUT`)
  CHRON_subgraph2tree=${L[2]}
fi


Radius=(`grep '^Ave. radius: ' $OUT`)


echo $VER ${Leaves[2]} ${Dissim[2]} ${CHRON_init[3]} ${CHRON_optim[4]} $CHRON_getBestChange $CHRON_tree2subgraph $CHRON_tree2subgraphDissim $CHRON_subgraph $CHRON_subgraph2tree ${Radius[2]}
