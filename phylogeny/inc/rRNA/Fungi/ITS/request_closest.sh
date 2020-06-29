#!/bin/bash
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "$0"
  exit 1
fi
NEW_OBJ=$1

#set -o xtrace 


INC=`dirname $0`
TMP=`mktemp`  
#echo $TMP  

blastn  -db $INC/seq.fa  -query $INC/../seq/$NEW_OBJ  -strand plus  -task blastn  -outfmt '6 sseqid nident' > $TMP.blastn
  # -evalue 1e-20
set +o errexit
grep -v $NEW_OBJ $TMP.blastn > $TMP.grep 
sort -k 2 -n -r $TMP.grep | cut -f 1 | head -100 > $TMP.head
set -o errexit
sort -u $TMP.head 


rm -fr $TMP*  



