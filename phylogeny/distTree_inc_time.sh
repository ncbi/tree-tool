#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: Distance tree data"
  echo "#2: Tree version"
  echo "#3: Chronometer::enabled (0/1)"
  exit 1
fi



Leaves=(`grep '^# Leaves: ' $1/hist/makeDistTree.$2`)
Dissim=(`grep '^# Dissimilarities: ' $1/hist/makeDistTree.$2`)
CHRON_topology=(`grep '^CHRON: Initial topology: ' $1/hist/makeDistTree.$2`)
CHRON_Leaves=(`grep '^CHRON: Optimizing new leaves: ' $1/hist/makeDistTree.$2`)


CHRON_getBestChange=0
CHRON_tree2subgraph=0
CHRON_tree2subgraphDissim=0
CHRON_subgraph=0
CHRON_subgraph2tree=0
if [ $3 == 1 ]; then
  L=(`grep '^CHRON: getBestChange: ' $1/hist/makeDistTree.$2`)
  CHRON_getBestChange=${L[2]}
  
  L=(`grep '^CHRON: tree2subgraph: ' $1/hist/makeDistTree.$2`)
  CHRON_tree2subgraph=${L[2]}
  
  L=(`grep '^CHRON: tree2subgraphDissim: ' $1/hist/makeDistTree.$2`)
  CHRON_tree2subgraphDissim=${L[2]}
  
  L=(`grep '^CHRON: subgraphOptimize: ' $1/hist/makeDistTree.$2`)
  CHRON_subgraph=${L[2]}
  
  L=(`grep '^CHRON: subgraph2tree: ' $1/hist/makeDistTree.$2`)
  CHRON_subgraph2tree=${L[2]}
fi


Radius=(`grep '^Ave. radius: ' $1/hist/makeDistTree.$2`)


#     #Leaves   #Dissim           topology           Leaves           getBestChange        tree2subgraph        tree2subgraphDissim        subgraph        subgraph2tree  Radius   
echo ${Leaves[2]} ${Dissim[2]} ${CHRON_topology[3]} ${CHRON_Leaves[4]} $CHRON_getBestChange $CHRON_tree2subgraph $CHRON_tree2subgraphDissim $CHRON_subgraph $CHRON_subgraph2tree ${Radius[2]}
