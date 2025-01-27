#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Print matching protein pairs"
  echo "#1: query protein sequences"
  echo "#2: subject protein sequences"
  echo "#3: min. fraction of identity"
  echo "#4: min. fraction of coverage"
  echo "#5: cores"
  exit 1
fi
QFASTA=$1
SFASTA=$2
IDENT=$3
COV=$4
CORES=$5


TMP=$( mktemp )
#comment $TMP 
#set -x


makeblastdb  -in $SFASTA   -dbtype prot  -blastdb_version 4  -logfile /dev/null  -out $TMP

blastp  -query $QFASTA  -db $TMP  -comp_based_stats 0  -evalue 1e-10  -seg no  -dbsize 10000000  -max_target_seqs 100000 -outfmt '6 qseqid sseqid qlen slen nident length qstart qend sstart send'  -num_threads $CORES > $TMP.blastp
#                                                                                                                                   1      2      3    4    5      6      7      8    9      10
awk '$5 >= $6 * '$IDENT' && ($8 - $7 + 1) >= $3 * '$COV' && ($10 - $9 + 1) >= $4 * '$COV'' $TMP.blastp | cut -f 1,2 


rm $TMP*
