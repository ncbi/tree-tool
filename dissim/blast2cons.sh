#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print the conservation dissimilarity of 2 DNA sequences"
  echo "#1: DNA 1 (FASTA)"
  echo "#2: DNA 2 (FASTA)"
  exit 1
fi
DNA1=$1
DNA2=$2


TMP=`mktemp`
#echo $TMP 


N1=`$THIS/../genetics/fasta2len $DNA1 -noprogress | cut -f 2`
N2=`$THIS/../genetics/fasta2len $DNA2 -noprogress | cut -f 2`

cp $DNA1 $TMP
makeblastdb  -in $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null
blastn  -db $TMP  -query $DNA2   -evalue 1e-100  -dust no  -task megablast  \
  -max_target_seqs 100000  -xdrop_gap 150  -xdrop_gap_final 150  -penalty -1  -gapopen 3  -gapextend 1  -dbsize 10000000  -searchsp 10000000 \
  -outfmt '6 qseqid sseqid qstart qend sstart send nident length' \
  | $THIS/blast2cons 0 $N2 $N1


rm $TMP*
