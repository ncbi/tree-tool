#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find closest virus"
  echo "#1: Virus.id"
  echo "#2: directory or ''"
  exit 1
fi
NEW_OBJ=$1
DIR="$2"


#set -o xtrace 


if [ $DIR ]; then
  error "Not implemented"
fi


TMP=`mktemp`  
#echo $TMP  


INC=`dirname $0`

blastn  -db $INC/seq.fa  -query $INC/../seq/$NEW_OBJ  -strand plus  -evalue 1e-20  -task blastn  -outfmt '6 sseqid nident' > $TMP.blastn
set +o errexit
grep -v $NEW_OBJ $TMP.blastn | sort -k 2 -n -r | cut -f 1 | head -100 > $TMP.head
set -o errexit
sort -u $TMP.head  


rm -fr $TMP*  



