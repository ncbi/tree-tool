#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Tight protein clustering"
  echo "#1: protein FASTA"
  echo "#2: cores"
  echo "#3: output file of representative sequences"
  exit 1
fi
F=$1
CORES=$2
OUT=$3


TMP=$( mktemp )
comment $TMP


$THIS/prot_match.sh $F $F 0.8 0.8 $CORES > $TMP.match
  # PAR
$THIS/fa2list.sh $TMP.fa | awk '{OFS="\t"; print $1, $1};' >> $TMP.match
$THIS/../connectPairs $TMP.match $TMP.pairs  -pairs -center
cut -f 2 $TMP.pairs > $OUT


rm $TMP*  
