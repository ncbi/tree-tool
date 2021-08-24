#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 3 ]; then
  echo "Contig coverage report"
  echo "#1: query FASTA file"
  echo "#2: subject FASTA file"
  echo "#3: dna_coverage mode: all, best, missed, combined"
  exit 1
fi
QUERY=$1
SUBJ=$2
MODE=$3


if [ ! -e $QUERY -o -d $QUERY ]; then
  error "Query file $QUERY does not exist"
fi

if [ ! -e $SUBJ -o -d $SUBJ ]; then
  error "Subject file $SUBJ does not exist"
fi

MT_MODE=""
N=`grep -c ">" $QUERY`
if [ $N == 0 ]; then
  error "$QUERY is not FASTA"
fi
if [ $N -gt 1 ]; then
  MT_MODE="-mt_mode 1"
fi


TMP=`mktemp`  
#echo $TMP > /dev/stderr
#set -x


echo "Running BLAST ..." > /dev/stderr
# DB
if [ -e $SUBJ.nhr ]; then
  DB=$SUBJ
else
  makeblastdb  -in $SUBJ  -out $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
  DB=$TMP
fi

HEADER="qseqid sseqid length nident qstart qend qlen sstart send slen stitle"
#       1      2      3      4      5      6    7    8      9    10   11
blastn  -query $QUERY  -db $DB  -dust no  -evalue 1e-100  -dbsize 10000000  -outfmt "6 $HEADER"  -num_threads 16  $MT_MODE | sort > $TMP.blastn
  # PAR

QUERY_NAME=`basename $QUERY`  
SUBJ_NAME=`basename $SUBJ`
$THIS/dna_coverage $TMP.blastn  -mode $MODE  -query "$QUERY_NAME"  -subject "$SUBJ_NAME"


rm $TMP*  
