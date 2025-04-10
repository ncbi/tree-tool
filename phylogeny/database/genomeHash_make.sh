#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Create GenomeHash.#3"
  echo "#1: genome/"
  echo "#2: #1 is large (0/1)"
  echo "#3: Suffix: HMM, PRT, CDS"
  exit 1
fi
GENOME=$1
LARGE=$2
SUF=$3


TMP=$( mktemp )


$THIS/../../check_file.sh $GENOME 0

# $TMP
if [ $LARGE -eq 1 ]; then
  $THIS/../../trav $GENOME "ls %d/%f" > $TMP
else
  ls $GENOME > $TMP
fi
wc -l $TMP

cp /dev/null GenomeHash.$SUF
$THIS/genomeHash_add.sh $GENOME $LARGE $TMP $SUF


rm $TMP*
