#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Output pairs of similar sequences"
  echo "#1; query   DNA FASTA"
  echo "#2; subject DNA blast daatbase"
  echo "#3: min. fraction of identity"
  echo "#4: min. fraction of coverage"
  echo "#5: blast num_threads"
  exit 1
fi
QUERY=$1
SUBJ=$2
IDENT=$3
COV=$4
THREADS=$5


blastn  -query $QUERY  -db $SUBJ  -dust no  -evalue 1e-10  -dbsize 10000000  -max_target_seqs 1000000 -outfmt '6 qseqid sseqid qlen slen nident length qstart qend sstart send'  -num_threads $THREADS \
  | awk '$5 >= $6 * '$IDENT' && ($8 - $7 + 1) >= $3 * '$COV' && ($10 > $9 && ($10 - $9 + 1) >= $4 * '$COV' || $10 < $9 && ($9 - $10 + 1) >= $4 * '$COV')' \
  | cut -f 1,2

