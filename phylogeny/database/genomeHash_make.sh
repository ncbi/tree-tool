#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Create GenomeHash.{CDS,PRT,HMM}"
  echo "#1: distance tree directory"
  echo "#2: #1/../genome is large (0/1)"
  echo "#3: Create Genomehash.CDS (0/1)"
  exit 1
fi
INC=$1
LARGE=$2
CDS=$3


TMP=$( mktemp )


GENOME=$INC/../genome
$THIS/../../check_file.sh $GENOME 0


# $TMP
if [ $LARGE -eq 1 ]; then
  $THIS/../../trav $GENOME "ls %d/%f" > $TMP
else
  ls $GENOME > $TMP
fi
wc -l $TMP

function run
{
  local SUF=$1
  H=""
  if [ $LARGE -eq 1 ]; then
    H="%h/"
  fi
  section "GenomeHash.$SUF"
  $THIS/../../trav $TMP "cat $GENOME/$H%f/%f.hash-$SUF | sed 's/^/%f\t/1'" > GenomeHash.$SUF
}

if [ $CDS -eq 1 ]; then
  run "CDS"
fi
run "PRT"
run "HMM"


rm $TMP*
