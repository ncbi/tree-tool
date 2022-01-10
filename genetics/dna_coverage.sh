#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 4 ]; then
  echo "Contig coverage report"
  echo "#1: query FASTA file"
  echo "#2: subject FASTA file"
  echo "#3: dna_coverage mode: all, best, missed, combine"
  echo "#4: force 1 line after header: 0/1"
  exit 1
fi
QUERY=$1
SUBJ=$2
MODE=$3
ONE_LINE=$4


#set -x


if [ ! -e $QUERY -o -d $QUERY ]; then
  error "Query file $QUERY does not exist"
fi

if [ ! -e $SUBJ -o -d $SUBJ ]; then
  error "Subject file $SUBJ does not exist"
fi

TMP=`mktemp`  
#echo $TMP > /dev/stderr


QUERY_NAME=`basename $QUERY`  
SUBJ_NAME=`basename $SUBJ`


# $TMP.blastn
N=`grep -c ">" $QUERY || true`
if [ $N -gt 0 -a -s $SUBJ ]; then
  MT_MODE=""
  if [ $N -gt 1 ]; then
    MT_MODE="-mt_mode 1"
  fi
  file $SUBJ > $TMP.file
  if grep ": gzip compressed data" $TMP.file &> /dev/null; then
    gunzip -c $SUBJ > $TMP.subj
    SUBJ=$TMP.subj
  fi
 #echo "Running BLAST ..." > /dev/stderr
  # DB
  if [ -e $SUBJ.nhr ]; then
    DB=$SUBJ
  else
    makeblastdb  -in $SUBJ  -out $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
    DB=$TMP
  fi
  HEADER="qseqid sseqid length nident qstart qend qlen sstart send slen stitle"
  #       1      2      3      4      5      6    7    8      9    10   11
  blastn  -query $QUERY  -db $DB  -dust no  -evalue 1e-100  -dbsize 10000000  -outfmt "6 $HEADER"  -num_threads 16  $MT_MODE 2> /dev/null | sort > $TMP.blastn
    # PAR
else
  touch $TMP.blastn
fi

FORCE=""
if [ $ONE_LINE -eq 1 ]; then
  FORCE="-force"
fi
PIDENT_MIN=""
if [ $MODE == "missed" ]; then
  PIDENT_MIN="-pident_min 90"  # PAR
fi
$THIS/dna_coverage $TMP.blastn  -mode $MODE  -query "$QUERY_NAME"  -subject "$SUBJ_NAME"  $PIDENT_MIN  $FORCE


rm $TMP*  
