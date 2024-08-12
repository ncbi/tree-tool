#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find 100 closest sequences in the tree"
  echo "#1: new sequence"
  echo "#2: directory of #1 or ''"  # ignored ??
  echo "#3: subset of sequence id's"  # ignored ??
  exit 1
fi
NEW_OBJ=$1

#set -o xtrace 


TMP=`mktemp`


INC=`dirname $0`
blastn  -db $INC/seq.fa  -query $INC/../seq/$NEW_OBJ  -strand plus  -evalue 1e-20  -task blastn  -outfmt '6 sseqid nident' | sort -k 2 -n -r | cut -f 1 > $TMP
head -100 $TMP | sort -u


rm $TMP*
