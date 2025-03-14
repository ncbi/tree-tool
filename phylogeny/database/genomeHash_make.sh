#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Create GenomeHash.{CDS,PRT,HMM}"
  echo "#1: distance tree directory"
  echo "#2: #1/../genome is large (0/1)"
  echo "#3: Suffix: HMM, PRT, CDS"
  exit 1
fi
INC=$1
LARGE=$2
SUF=$3



TMP=$( mktemp )


GENOME=$INC/../genome
$THIS/../../check_file.sh $GENOME 0

# $TMP, H
if [ $LARGE -eq 1 ]; then
  $THIS/../../trav $GENOME "ls %d/%f" > $TMP
  H="%h/"
else
  ls $GENOME > $TMP
  H=""
fi
wc -l $TMP

$THIS/../../trav $TMP "cat $GENOME/$H%f/%f.hash-$SUF | sed 's/^/%f\t/1'" > GenomeHash.$SUF



rm $TMP*
