#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Append GenomeHash.#4"
  echo "#1: genome/"
  echo "#2: #1 is large (0/1)"
  echo "#3: list of genomes to add"
  echo "#4: Suffix: HMM, PRT, CDS"
  exit 1
fi
GENOME=$1
LARGE=$2
ADD=$3
SUF=$4


$THIS/../../check_file.sh $GENOME 0

OUT=GenomeHash.$SUF
if [ ! -e $OUT ]; then
  $THIS/genomeHash_make.sh $GENOME $LARGE $SUF
fi

H=""
if [ $LARGE -eq 1 ]; then
  H="%h/"
fi

$THIS/../../trav $ADD "cat $GENOME/$H%f/%f.hash-$SUF | sed 's/^/%f\t/1'" >> $OUT

