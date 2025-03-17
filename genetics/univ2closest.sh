#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
GOAL=100
if [ $# -ne 3 ]; then
  echo "Return: top $GOAL BLASTP hits"
  echo "#1; query sequences of universal proteins"
  echo "#2: subject protein BLAST database of universal proteins"
  echo "#3: restrict subject id's to this list | ''"
  exit 1
fi
QUERY=$1
DB=$2
TARGET=$3  


TMP=$( mktemp )


blastp  -query $QUERY  -db $DB  -task blastp-fast  -threshold 100  -window_size 15  -comp_based_stats 0  -seg no  -evalue 1e-10  -num_threads 5  -outfmt '6 sseqid slen nident' | tr '-' '\t' | awk '{OFS="\t"; print $1, $3 - $4};' | sort -k2,2n | cut -f 1 | uniq > $TMP
if [ ! -s $TMP ]; then
  error "Empty BLASTP output"
fi

# $TMP.target: ordered list with duplications
if [ "$TARGET" ]; then
  $THIS/../check_file.sh $TARGET 1
  grep -xf $TARGET $TMP > $TMP.target || true
else
  mv $TMP $TMP.target
fi

$THIS/../multilist2subset $TMP.target $GOAL


rm $TMP*
