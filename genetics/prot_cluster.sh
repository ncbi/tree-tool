#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Cluster protein sequences, if #1 = #2 then print representatives else print #1-#2 pairs"
  echo "#1: query directory with protein sequence"
  echo "#2: subject directory with protein sequence"
  echo "#3: min. fraction of identity"
  echo "#4: min. fraction of coverage"
  echo "#5: output file with pairs of FASTA identfiers"
  exit 1
fi
QDIR=$1
SDIR=$2
IDENT=$3
COV=$4
PAIR=$5


TMP=$( mktemp )
#comment $TMP 
#set -x


$THIS/../trav $SDIR "cat %d/%f" > $TMP
makeblastdb  -in $TMP   -dbtype prot  -blastdb_version 4  -logfile /dev/null

# $PAIR
if [ $QDIR == $SDIR ]; then
  $THIS/fa2list.sh $TMP | awk '{OFS="\t"; print $1, $1};' > $PAIR
fi
$THIS/../trav -step 1 $QDIR "blastp  -query %d/%f  -db $TMP  -comp_based_stats 0  -evalue 1e-10  -seg no  -dbsize 10000000  -max_target_seqs 100000 -outfmt '6 qseqid sseqid qlen slen nident length qstart qend sstart send'  -num_threads 16" > $TMP.blastp
#                                                                                                                                                              1      2      3    4    5      6      7      8    9      10
awk '$5 >= $6 * '$IDENT' && ($8 - $7 + 1) >= $3 * '$COV' && ($10 - $9 + 1) >= $4 * '$COV'' $TMP.blastp | cut -f 1,2 >> $PAIR

if [ $QDIR == $SDIR ]; then
  $THIS/../connectPairs $PAIR $TMP.out  -pairs 
  cut -f 2 $TMP.out | sort -u
else
  cat $PAIR
fi


rm $TMP*
