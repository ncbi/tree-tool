#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print non-standard start codons (not 'M')"
  echo "#1: Peptide FASTA-file"
  exit 1
fi

cat $1 | tr '\n' ' ' | tr '>' '\n' | grep -v '^$' | sed 's/^\([0-9]\+\) \(.\).*$/\1 \2/1' | grep -v ' M'
