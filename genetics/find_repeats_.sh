#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find repeats for a fragment"
  echo "#1: DNA fragment"
  echo "#2: DNA BLAST DB"
  echo "#3: hanging tail length"
  echo "#4: output durectory"
  exit 1
fi
FRAG=$1
DB=$2
HANG=$3
OUT=$4


TMP=$( mktemp )


blastn  -query $FRAG  -db $DB  -evalue 1e-20  -dust no  -dbsize 10000000  -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen' > $TMP.blastn
  #                                                                                  1      2      3      4      5      6    7    8      9    10
awk -v H=$HANG '$4 / $3 >= 0.9 && $5 <= H && $7 - $6 <= H' $TMP.blastn > $TMP.blastn1
N=$( < $TMP.blastn1  wc -l )
if [ $N -gt 2 ]; then  # PAR
  cut -f 2,8,9 $TMP.blastn1 > $OUT/$( basename $FRAG )
fi


rm $TMP*
