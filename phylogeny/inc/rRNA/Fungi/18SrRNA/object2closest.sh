#!/bin/bash --noprofile
source $PANFS/code/cpp/bash_common.sh
if [ $# -ne 1 ]; then
  echo "$0"
  exit 1
fi
NEW_OBJ=$1

#set -o xtrace 


TMP=`mktemp`


INC=`dirname $0`
blastn  -db $INC/seq.fa  -query $INC/../seq/$NEW_OBJ  -strand plus  -evalue 1e-20  -task blastn  -outfmt '6 sseqid nident' | sort -k 2 -n -r | cut -f 1 > $TMP
head -100 $TMP | sort -u


rm $TMP*
