#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Output: #3/<name of #2>"
  echo "#1: Protein sequences"
  echo "#2: DNA sequences"
  echo "#3: Output directory"
  exit 1
fi
PROT=$1
NUCL=$2
OUT=$3


TMP=$( mktemp )
#echo $TMP


NAME=$( basename $PROT )

cp $PROT $TMP.fa
makeblastdb  -in $TMP.fa  -dbtype prot  -logfile /dev/null

blastx  -db $TMP.fa  -query $NUCL  -show_gis  -word_size 6  -evalue 1e-10  -query_gencode 11  -seg no  -comp_based_stats 0  -xdrop_gap_final 1000  -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen' > $OUT/$NAME
#                                                                                                                                                             1      2      3      4      5      6    7     8     9    10


rm -f $TMP*
