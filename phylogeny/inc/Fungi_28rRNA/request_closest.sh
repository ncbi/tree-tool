#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "$0"
  exit 1
fi
NEW_OBJ=$1

#set -o xtrace 


BASE_DIR=`dirname $0`
TMP=`mktemp`  

blastn  -db $BASE_DIR/seq.fa  -query /home/brovervv/panfs/marker/Fungi/28SrRNA/seq/$NEW_OBJ  -strand plus  -evalue 1e-20  -outfmt '6 sseqid' > $TMP.blastn
set +o errexit
grep -v $NEW_OBJ $TMP.blastn | sort -u > $TMP.grep
set -o errexit
head -100 $TMP.grep | sed 's/$/ '$NEW_OBJ'/1'


rm -fr $TMP*



