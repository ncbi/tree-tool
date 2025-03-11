#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Splice introns out"
  echo "#1: input protein multi-FASTA file of similar proteins"
  echo "#2: output multi-FASTA file"
  echo "#3: cores"
  exit 1
fi
F=$1
OUT=$2
CORES=$3


TMP=$( mktemp )
comment $TMP 


mkdir $TMP.seq
$THIS/splitFasta $F $TMP.seq -aa -qc 

mkdir $TMP.ex
$THIS/../trav $TMP.seq "$THIS/blastp_seq2exons.sh %d/%f $F  1> $TMP.ex/%f  2> /dev/null"  -threads $CORES  -step 1
$THIS/../trav $TMP.ex "cat %d/%f" > $OUT


rm -r $TMP*  
