#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Set #1.prot-univ"
  echo "#1: Genome.id"
  echo "#2: large (0/1)"
  echo "#3: suffix: HMM|tblastn"
  exit 1
fi
G=$1
LARGE=$2
SUF=$3


# $F
H=""
if [ $LARGE -eq 1 ]; then
  H=$( $THIS/../../file2hash $G )
fi
F=genome/$H/$G/$G.prot-univ

$THIS/../../check_file.sh $F.$SUF 1

rm -f $F
ln -s $PWD/$F.$SUF $F


