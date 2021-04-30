#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# != 1 ]; then
  echo "#1: Genome.id"
  exit 1
fi
G=$1


H=`$THIS/../../file2hash $G`
N=`cat genome/$H/$G/$G.hash-CDS | wc -l`
echo $G $N
