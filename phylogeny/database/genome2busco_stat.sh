#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Input: genome/<hash>/#1/#1.busco10.tsv" 
  echo "Print: <Genome.id> <#complete_genes> <#duplicate_genes> <#total_genes> <duplicate_genes_share>"
  echo "#1: Genome.id"
  exit 1
fi
G=$1


TMP=$( mktemp )


H=$( $THIS/../../file2hash $G )

grep -v '^#' genome/$H/$G/$G.busco10.tsv | cut -f 1,2 > $TMP
C=$( awk -F '\t' '$2 == "Complete"'   $TMP | cut -f 1 | sort -u | wc -l )
D=$( awk -F '\t' '$2 == "Duplicated"' $TMP | cut -f 1 | sort -u | wc -l )
F=$( awk -F '\t' '$2 == "Fragmented"' $TMP | cut -f 1 | sort -u | wc -l )
echo -e "$G\t$C\t$D\t$F"


rm $TMP
