#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find frame shifts"
  echo "#1: DNA (FASTA)"
  echo "#2: reference proteins (FASTA)"
  echo "#3: protens end with stop codons (0/1)"
  echo "#4: gencode"
  exit 1
fi
DNA=$1
PROT=$2
STOP_CODON=$3
GC=$4


TMP=$( mktemp )
comment $TMP
#set -x


makeblastdb  -in $DNA  -dbtype nucl  -out $TMP  -blastdb_version 4  -logfile /dev/null
# PAR
tblastn  -db $TMP  -query $PROT  -db_gencode $GC  -seg no  -comp_based_stats 0  -max_target_seqs 10000  -word_size 3  -threshold 21  -dbsize 10000  -evalue 1  -num_threads 15  -outfmt '6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq' | sort > $TMP.blast
STOP_CODON_PAR=""
if [ $STOP_CODON == 1 ]; then
  STOP_CODON_PAR="-stop_codon"
fi
$THIS/tblastn2disruption $TMP.blast -noprogress $STOP_CODON_PAR -bacteria -qc > $TMP.fs
$THIS/disruption2genesymbol $DNA $PROT $TMP.fs  -gencode $GC  -noprogress -qc


rm $TMP*
