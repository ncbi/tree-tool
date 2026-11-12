#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: BLAST output: sseqid bitscore"
  exit 1
fi
B=$1  


TMP=$( mktemp )


TOP=$( sort -k2nr $B | head -1 | cut -f 2 )
THRESHOLD=$( echo "$TOP * 0.9 / 1" | bc )  # PAR
awk -F '\t' '$2 >= '$THRESHOLD $B > $TMP
N=$( < $TMP  wc -l )
sed 's/$/\t'$N'/1' $TMP


rm $TMP*
