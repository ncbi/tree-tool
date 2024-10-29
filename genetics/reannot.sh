#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Find ORFs by tblastn"
  echo "#1: DNA multi-FASTA file"
  echo "#2: reference proteins"
  echo "#3: gencode"
  echo "#4: output protein file"
  echo "#5: log-file"
  exit 1
fi
DNA=$1
REF=$2
GENCODE=$3
PROT=$4
LOG=$5


TMP=$( mktemp )
#echo $TMP 


makeblastdb  -in $DNA  -dbtype nucl  -blastdb_version 4  -out $TMP  -logfile /dev/null 
# PAR
tblastn  -db $TMP  -query $REF  -show_gis  -word_size 3  -evalue 1e-4  -db_gencode $GENCODE  -seg no  -comp_based_stats 0  -outfmt '6 qseqid sseqid length positive qstart qend sstart send slen sseq' | awk '$4/$3 >= 0.85' > $TMP.tblastn
$THIS/tblastn2orfs  -noprogress  $DNA $TMP.tblastn $GENCODE  -log $LOG > $PROT


rm -f $LOG
rm $TMP*
