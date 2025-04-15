#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find frame shifts"
  echo "#1: DNA (FASTA)"
  echo "#2: reference proteins (FASTA)"
  echo "#3: gencode"
  exit 1
fi
DNA=$1
PROT=$2
GC=$3


TMP=$( mktemp )
#comment $TMP
#set -x


makeblastdb  -in $DNA  -dbtype nucl  -out $TMP  -blastdb_version 4  -logfile /dev/null
# PAR
tblastn  -db $TMP  -query $PROT  -show_gis  -db_gencode $GC  -seg no  -comp_based_stats 0  -max_target_seqs 10000  -word_size 3  -threshold 21  -num_threads 15  -outfmt '6 qseqid sseqid qstart qend sstart send qseq sseq' | sort > $TMP.blast
$THIS/tblastn2frameshift $TMP.blast


rm $TMP*
