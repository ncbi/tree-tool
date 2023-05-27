#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Return: top 100 BLASTP hits"
  echo "#1; query sequences of universal proteins"
  echo "#2: subject protein BLAST database of universal proteins"
  exit 1
fi
QUERY=$1
DB=$2


TMP=`mktemp`


blastp  -db $DB  -query $QUERY  -task blastp-fast  -threshold 100  -window_size 15  -comp_based_stats 0  -seg no  -evalue 1e-10  -num_threads 5  -outfmt '6 sseqid slen nident' | tr '-' '\t' | awk '{OFS="\t"; print $1, $3 - $4};' | sort -k2,2n | cut -f 1 | uniq > $TMP
head -100 $TMP | sort -u


rm $TMP*
