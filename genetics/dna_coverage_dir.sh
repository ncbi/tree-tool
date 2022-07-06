#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 6 ]; then
  echo "Contig coverage report for a directory of queries"
  echo "#1: diretory with query FASTA files"
  echo "#2: subject FASTA file (can be gzip'ped or converted into a BLAST database)"
  echo "#3: dna_coverage mode: all, combine"
  echo "#4: min. identity percent"
  echo "#5: min. alignment length"
  echo "#6: output report"
  exit 1
fi
QUERY_DIR=$1
SUBJ=$2
MODE=$3
PIDENT_MIN=$4
ALIGN_MIN=$5
OUT=$6



TMP=`mktemp`


$THIS/../trav $QUERY_DIR  -step 1  "dna_coverage.sh %d/%f $SUBJ $MODE $PIDENT_MIN $ALIGN_MIN '' 0" > $TMP
$THIS/../tsv/tsv_clean.sh $TMP > $OUT


rm $TMP*
