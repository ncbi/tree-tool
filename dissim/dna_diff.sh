#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print the difference of 2 DNA sequences by semiglobal alignment"
  echo "#1: DNA 1 (FASTA)"
  echo "#2: DNA 2 (FASTA)"
  echo "#3: min. alignment length"
  exit 1
fi
DNA1=$1
DNA2=$2
LEN_MIN=$3


TMP=`mktemp`
#echo $TMP 


cp $DNA1 $TMP
makeblastdb  -in $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null
set +o pipefail
blastn  -db $TMP  -query $DNA2   -evalue 1e-100  -dust no  -task blastn  -outfmt '6 length qstart qend qlen sstart send slen qseq sseq' | head -1 | awk '$1 > '$LEN_MIN' && ($2 == 1 || $5 == 1) && ($3 == $4 || $6 == $7)' | awk '{printf "%s\n%s\n", $8, $9};' > $TMP.pair
#                                                                                   1      2      3    4    5      6    7    8    9    
if [ -s $TMP.pair ]; then
  $THIS/dna_diff $TMP.pair 
else
  echo "inf"
fi


rm $TMP*
