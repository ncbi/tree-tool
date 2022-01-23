#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 5 ]; then
  echo "Contig coverage report"
  echo "#1: query FASTA file"
  echo "#2: subject FASTA file"
  echo "#3: dna_coverage mode: all, best, missed, combine"
  echo "#4: sequence identifier to be removed from #2 or ''"
  echo "#5: force 1 line after header: 0/1"
  exit 1
fi
QUERY=$1
SUBJ=$2
MODE=$3
TO_REMOVE="$4"
ONE_LINE=$5


#set -x


if [ ! -e $QUERY -o -d $QUERY ]; then
  error "Query file $QUERY does not exist"
fi


TMP=`mktemp`  
#echo $TMP > /dev/stderr


QUERY_NAME=`basename $QUERY`  
SUBJ_NAME=`basename $SUBJ`


DB=$SUBJ
if [ ! -e $DB.nhr ]; then
  if [ ! -e $SUBJ -o -d $SUBJ ]; then
    rm $TMP*
    error "Subject file $SUBJ does not exist"
  fi
  file $SUBJ > $TMP.file
  if grep ": gzip compressed data" $TMP.file &> /dev/null; then
    gunzip -c $SUBJ > $TMP.subj
    SUBJ=$TMP.subj
    DB=$SUBJ
  fi
  if [ -s $SUBJ ]; then
    makeblastdb  -in $SUBJ  -out $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
    DB=$TMP
  fi
fi


# $TMP.blastn
N=`grep -c ">" $QUERY || true`
if [ $N -gt 0 -a -e $DB.nhr ]; then
  MT_MODE=""
  if [ $N -gt 1 ]; then
    MT_MODE="-mt_mode 1"
  fi
 #echo "Running BLAST ..." > /dev/stderr
  HEADER="qseqid sseqid length nident qstart qend qlen sstart send slen stitle"
  #       1      2      3      4      5      6    7    8      9    10   11
  blastn  -query $QUERY  -db $DB  -dust no  -evalue 1e-10  -dbsize 10000000  -outfmt "6 $HEADER"  -num_threads 16  $MT_MODE 2> /dev/null | sort > $TMP.blastn
    # PAR
  if [ $TO_REMOVE ]; then
    awk '$2 != "'$TO_REMOVE'"' $TMP.blastn > $TMP.blastn1
    mv $TMP.blastn1 $TMP.blastn
  fi
else
  cp /dev/null $TMP.blastn
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
