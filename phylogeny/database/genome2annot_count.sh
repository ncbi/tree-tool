#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# != 2 ]; then
  echo "#1: Genome.id"
  echo "#2: CDS|PRT|HMM"
  exit 1
fi
G=$1
SUF=$2


H=$( $THIS/../../file2hash $G )
N=$( < genome/$H/$G/$G.hash-$SUF  wc -l )
echo $G $N
