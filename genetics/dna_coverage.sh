#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# != 8 ]; then
  echo "Contig coverage report"
  echo "#1: query DNA FASTA file"
  echo "#2: subject DNA FASTA file (can be gzip'ped or converted into a BLAST database)"
  echo "#3: dna_coverage mode: all, combine"
  echo "#4: min. identity percent of HSP"
  echo "#5: min. alignment length of HSP"
  echo "#6: sequence identifier to be removed from #2 or ''"
  echo "#7: force 1 line after header: 0/1"
  echo "#8: include all query sequences: 0/1"
  exit 1
fi
QUERY=$1
SUBJ=$2
MODE=$3
PIDENT_MIN=$4
ALIGN_MIN=$5
TO_REMOVE="$6"
ONE_LINE=$7
ALL_QUERY=$8


#set -x


if [ ! -e $QUERY ] || [ -d $QUERY ]; then
  error "Query file $QUERY does not exist"
fi


TMP=$( mktemp )
#comment $TMP 
#set -x


QUERY_NAME=$( basename $QUERY )
SUBJ_NAME=$( basename $SUBJ )


DB=$SUBJ
DB_EXISTS=0
if [ -e $DB.nhr ] || [ -e $DB.00.nhr ]; then
  DB_EXISTS=1
fi
if [ $DB_EXISTS == 0 ]; then
  if [ ! -e $SUBJ ] || [ -d $SUBJ ]; then
    rm $TMP*
    error "Subject file $SUBJ does not exist"
  fi
  file $SUBJ > $TMP.file
  if grep ": gzip compressed data" $TMP.file &> /dev/null; then
    gunzip -c $SUBJ > $TMP.subj
    SUBJ=$TMP.subj
  fi
  if [ -s $SUBJ ]; then
    makeblastdb  -in $SUBJ  -out $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
    DB=$TMP
    DB_EXISTS=1
  fi
fi

# $TMP.blastn
N=$( grep -c ">" $QUERY || true )
if [ $N -gt 0 ] && [ $DB_EXISTS == 1 ]; then
  MT_MODE=""
  if [ $N -gt 1 ]; then
    MT_MODE="-mt_mode 1"
  fi
  HEADER="qseqid sseqid length nident qstart qend qlen sstart send slen stitle"
  #       1      2      3      4      5      6    7    8      9    10   11
  blastn  -query $QUERY  -db $DB  -dust no  -max_target_seqs 1000000  -outfmt "6 $HEADER"  -num_threads 16  2> /dev/null | sort -k1,1 -k4,4nr > $TMP.blastn
    # $MT_MODE ??
    # PAR
  if [ "$TO_REMOVE" ]; then
    awk '$2 != "'$TO_REMOVE'"' $TMP.blastn > $TMP.blastn1
    mv $TMP.blastn1 $TMP.blastn
  fi
else
  > $TMP.blastn
fi

FORCE=""
if [ $ONE_LINE -eq 1 ]; then
  FORCE="-force"
fi
$THIS/dna_coverage $TMP.blastn  -mode $MODE  -query "$QUERY_NAME"  -subject "$SUBJ_NAME"  -pident_min $PIDENT_MIN  -align_min $ALIGN_MIN  $FORCE > $TMP.cov

if [ $ALL_QUERY == 1 ]; then
  echo -e "#query\tqseqid\tqlen" > $TMP.len
  $THIS/fasta2len $QUERY -noprogress | sed 's/^/'$QUERY_NAME'\t/1' >> $TMP.len
  $THIS/../tsv/tsv_expand.sh $TMP.len $TMP.cov '-left' &> /dev/null
  cat $TMP.len
else
  cat $TMP.cov
fi


rm $TMP*  
