#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Return: top 100 BLASTN hits in top strand"
  echo "#1; query DNA sequence"
  echo "#2: subject DNA BLAST database"
  echo "#3: subset of sequence id's or ''"
  exit 1
fi
QUERY=$1
DB=$2
SUBSET=$3


TMP=$( mktemp )
#comment $TMP


blastn  -db $DB  -query $QUERY  -strand plus  -task blastn  -num_threads 5  -outfmt '6 sseqid nident' | sort -u > $TMP

# $TMP -> $TMP.inter
if [ "$SUBSET" ]; then
  sort -cu $SUBSET
  join  -1 1  -2 1  $TMP $SUBSET | tr ' ' '\t' > $TMP.inter
else
  mv $TMP $TMP.inter 
fi

sort -k2nr $TMP.inter | cut  -f 1  > $TMP.sorted
head -100 $TMP.sorted 


rm $TMP*
