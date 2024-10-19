#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print the CDSs of FASTA file with wrong start/stop codons"
  echo "#1: CDSs FASTA file"
  echo "#2: Output file"
  exit 1
fi
FASTA=$1
OUT=$2


TMP=$( mktemp )
comment $TMP


mkdir $TMP.seq
$THIS/splitFasta $FASTA $TMP.seq

$THIS/../trav -step 1  -errors /dev/stdout  $TMP.seq  "CDS_check.sh %d/%f" > $OUT


rm -r $TMP*
