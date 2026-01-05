#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Output: #3/<name of #2>"
  echo "#1: protein FASTA"
  echo "#2: DNA FASTA (can be gzip'ed)"
  echo "#3: output is tabular (0/1)"
  echo "#4: output file"
  exit 1
fi
PROT=$1
ASM=$2
TAB=$3
OUT=$4


TMP=$( mktemp )
#comment $TMP


EXT=$( echo $ASM | tr '.' '\n' | tail -1 )
if [ $EXT == "gz" ]; then
  cp $ASM $TMP.gz
  rm $TMP
  gunzip $TMP.gz
else
  cp $ASM $TMP
fi
makeblastdb  -in $TMP  -dbtype nucl  -logfile /dev/null

OUTFMT=""
if [ $TAB == 1 ]; then
  OUTFMT="-outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen'"
  #                  1      2      3      4       5     6    7    8      9    10
fi

tblastn  -db $TMP  -query $PROT  -show_gis  -word_size 6  -evalue 1e-10  -db_gencode 11  -seg no  -comp_based_stats 0  -xdrop_gap_final 1000  $OUTFMT  > $OUT


rm -f $TMP*
