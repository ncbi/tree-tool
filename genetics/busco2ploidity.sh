#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print: file #complete_genes #duplicate_genes #total_genes duplicate_genes_share"
  echo "#1: BUSCO .tsv-file"
  exit 1
fi
IN=$1


TMP=$( mktemp )


grep -v '^#' $IN | cut -f 1,2 > $TMP
C=$( awk -F '\t' '$2 == "Complete"'   $TMP | cut -f 1 | sort -u | wc -l )
D=$( awk -F '\t' '$2 == "Duplicated"' $TMP | cut -f 1 | sort -u | wc -l )
T=$(( C + D ))
if [ $T == 0 ]; then
  F="?"
else
  F=$( echo "$D / $T" | bc -l )
fi
echo -e "$IN\t$C\t$D\t$T\t$F"


rm $TMP
