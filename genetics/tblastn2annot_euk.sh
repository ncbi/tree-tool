#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find proteins with introns by tblastn"
  echo "#1: DNA multi-FASTA file"
  echo "#2: reference proteins"
  echo "#3: output protein file"
  exit 1
fi
DNA=$1
REF=$2
PROT=$3


TMP=$( mktemp )
comment $TMP  
#set -x  

makeblastdb  -in $DNA  -out $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
tblastn  -query $REF  -db $TMP  -db_gencode 1  -evalue 1000  -word_size 2  -outfmt '6 qseqid sseqid qstart qend sstart send qseq sseq' -num_threads 4 > $TMP.tblastn
  # -task tblastn-fast  -word_size 5  -evalue 1e-4: weak HSPs are needed!
  # -num_threads requires -db
  # -seg no  -comp_based_stats 0 
$THIS/tblastn2annot_euk  $TMP.tblastn  -qc  -log $TMP.log > $PROT


rm $TMP*  
