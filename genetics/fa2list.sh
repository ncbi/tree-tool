#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print the list of ids of a FASTA file"
  echo "#1: FASTA-file"
  exit 1
fi
F=$1


grep '^>' $F | sed 's/^>//1'| sed 's/ .*$//1'
