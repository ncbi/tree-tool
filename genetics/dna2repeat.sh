#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 1 ]; then
  echo "Check if a DNA is repeated"
  echo "#1: DNA FASTA file"
  exit 1
fi
DNA=$1


TMP=`mktemp`  
#echo $TMP  


cp $DNA $TMP
makeblastdb  -in $TMP   -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
blastn  -db $TMP  -query $TMP  -outfmt '6 qstart qend sstrand sstart send length nident' | awk '$1 != $4 || $2 != $5'
#                                         1      2    3       4      5    6      7


rm -f $TMP*  
