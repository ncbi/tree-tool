#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Return: top #3 BLASTN hits in top strand"
  echo "#1; Query DNA sequence"
  echo "#2: Subject DNA database"
  echo "#3: Number of top hits"
  exit 1
fi
QUERY=$1
DB=$2
NUM=$3



TMP=`mktemp`  
#echo $TMP  


blastn  -db $DB  -query $QUERY  -strand plus  -task blastn  -outfmt '6 sseqid nident' > $TMP.blastn
NAME=`head -1 $QUERY | sed 's/^>//1' | cut -f 1 -d ' '`
set +o errexit
grep -v $NAME $TMP.blastn > $TMP.grep 
sort -k 2 -n -r $TMP.grep | cut -f 1 | head -$NUM > $TMP.head
set -o errexit
sort -u $TMP.head 


rm $TMP*  



