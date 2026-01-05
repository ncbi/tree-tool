#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print start codons"
  echo "#1: FASTA-file"
  echo "#2: DNA - 0; peptide - 1"
  exit 1
fi
F=$1
PEP=$2

DOTS="..."
if [ $PEP == 1 ]; then
  DOTS="."
fi
sed 's/ .*$//1' $F | tr '\n' ' ' | tr '>' '\n' | grep -v '^$' | sed 's/^\([^ ]\+\) \('$DOTS'\).*$/\1\t\2/1' 
