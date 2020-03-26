#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "$0"
  exit 1
fi
NEW_OBJ=$1

#set -o xtrace 


INC=`dirname $0`
TMP=`mktemp`  
#echo $TMP  

blastn  -db $INC/seq.fa  -query /home/brovervv/panfs/marker/bacteria/seq/$NEW_OBJ  -strand plus  -evalue 1e-20  -outfmt '6 sseqid nident' > $TMP.blastn
set +o errexit
grep -v $NEW_OBJ $TMP.blastn > $TMP.grep 
set -o errexit
sort -k 2 -n -r $TMP.grep | cut -f 1 | head -100 | sort -u | sed 's/$/ '$NEW_OBJ'/1' 


rm -fr $TMP*  



