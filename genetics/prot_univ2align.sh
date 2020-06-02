#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: Genome.id"
  exit 1
fi
OBJ=$1


rm -rf genome/$1/seq  # ??

mkdir genome/$1/seq
$THIS/splitFastaProt -log log/$1  genome/$1/$1.prot-univ 25 genome/$1/seq
$THIS/../trav genome/$1/seq -log log/$1 "hmmalign --amino --informat FASTA --outformat A2M hmm/%f.HMM %d/%f" | sed '/^[^>]/ s/[a-z]//g' | sed '/^[[:space:]]*$/d' > genome/$1/$1.hmm-align
rm -r genome/$1/seq

rm -f log/$1
