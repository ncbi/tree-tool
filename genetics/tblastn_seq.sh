#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Output: #3/<name of #2>"
  echo "#1: protein FASTA"
  echo "#2: DNA FASTA"
  echo "#3: output file"
  exit 1
fi
PROT=$1
ASM=$2
OUT=$3


TMP=`mktemp`
#echo $TMP


NAME=`basename $ASM`

cp $ASM $TMP.fa
makeblastdb  -in $TMP.fa  -dbtype nucl  -logfile /dev/null

tblastn  -db $TMP.fa  -query $PROT  -show_gis  -word_size 6  -evalue 1e-10  -db_gencode 11  -seg no  -comp_based_stats 0  -xdrop_gap_final 1000  -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen' > $OUT
#                                                                                                                                                           1      2      3      4       5     6    7    8      9    10


rm -f $TMP*
