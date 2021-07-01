#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 3 ]; then
  echo "Contig coverage report"
  echo "#1: query QUERY FASTA file"
  echo "#2: subject QUERY FASTA file"
  echo "#3: report non-covered query segments (0/1)"
  exit 1
fi
QUERY=$1
SUBJ=$2
MISSED=$3


# Threashold
T=90.0  


TMP=`mktemp`  
#echo $TMP > /dev/stderr
#set -x


HEADER="qseqid sseqid length nident qstart qend qlen sstart send slen stitle"
#       1      2      3      4      5      6    7    8      9    10   11

echo "Running BLAST ..." > /dev/stderr
makeblastdb  -in $SUBJ  -out $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
blastn  -query $QUERY  -db $TMP  -dust no  -evalue 1e-100  -dbsize 10000000  -outfmt "6 $HEADER"  -num_threads 16  -mt_mode 1 | sort > $TMP.blastn


if [ $MISSED == 1 ]; then
  dna_coverage $TMP.blastn  -mode missed
else
  dna_coverage $TMP.blastn | awk -F '\t' '$10 > '$T' && $16 > '$T
fi


rm $TMP*  
