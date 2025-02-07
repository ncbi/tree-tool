#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print #1 with introns spliced out"
  echo "#1: query protein sequence"
  echo "#2: subject protein database"
  exit 1
fi
Q=$1
DB=$2


TMP=$( mktemp )
#comment $TMP 


blastp  -query $Q  -subject $DB  -outfmt '6 qseqid sseqid qstart qend qlen qseq sseq'  -max_target_seqs 100000  > $TMP.blastp
$THIS/blastp2exons $TMP.blastp  


rm $TMP*  
