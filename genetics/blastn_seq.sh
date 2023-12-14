#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Output: #3/<name of #2>"
  echo "#1: query DNA sequences"
  echo "#2: subject DNA sequences"
  exit 1
fi
QUERY=$1
SUBJ=$2


TMP=`mktemp`
#echo $TMP


NAME=`basename $QUERY`

# $TMP.fa
GZ=`echo $SUBJ | tr '.' '\n' | tail -1`
if [ $GZ == "gz" ]; then
  gunzip -c $SUBJ > $TMP.fa
else
  cp $SUBJ $TMP.fa
fi

makeblastdb  -in $TMP.fa  -dbtype nucl  -logfile /dev/null

blastn  -db $TMP.fa  -query $QUERY  -show_gis  -evalue 1e-20  -dust no  -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen'
#                                                                                  1      2      3      4      5      6    7    8      9    10


rm -f $TMP*
