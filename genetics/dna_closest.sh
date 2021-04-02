#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Return: top 100 BLASTN hits in top strand"
  echo "#1; Query DNA sequence"
  echo "#2: Subject DNA database"
  exit 1
fi
QUERY=$1
DB=$2


blastn  -db $DB  -query $QUERY  -strand plus  -task blastn  -num_threads 5  -outfmt '6 sseqid nident' | sort -k 2 -n -r | cut -f 1 | head -100 | sort -u
