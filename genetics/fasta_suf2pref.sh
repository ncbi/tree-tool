#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print FASTA changing identifiers: PREF-SUF -> SUF-PREF, where SUF contains no '-'"
  echo "#1: input FASTA"
  exit 1
fi
F=$1


sed 's/^>\(.\+\)-\([^-]\+\)$/>\2-\1/1' $F