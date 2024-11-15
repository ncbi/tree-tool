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
#comment $TMP  
#set -x  


# PAR
tblastn  -query $REF  -subject $DNA  -db_gencode 1  -seg no  -comp_based_stats 0  -outfmt '6 qseqid sseqid qstart qend sstart send qseq sseq' > $TMP.tblastn
  # -task tblastn-fast  -word_size 5  -evalue 1e-4: weak HSPs are needed!
$THIS/tblastn2annot_euk  $TMP.tblastn  -qc  -log $TMP.log > $PROT


rm $TMP*  
