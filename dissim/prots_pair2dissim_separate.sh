#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo 'Input: genome/*/*/*.prot-univ, hmm-univ.list'
  echo "Output: #4/#1, #3.log/#1"
  echo "#1: file in #3 with pairs of genomes"
  echo "#2: 1 - BLOSUM62, 0 - PAM30"
  echo "#3: input directory"
  echo "#4: output directory"
  exit 1
fi
FILE=$1
BLOSUM62_ARG=$2
IN_DIR=$3
OUT_DIR=$4


GENOME=genome  

LOG=$IN_DIR.log/$FILE

BLOSUM62=""
if [ $BLOSUM62_ARG == 1 ]; then
  BLOSUM62="-blosum62"
fi


TMP=`mktemp`
#echo $TMP


rm -f $OUT_DIR/$FILE

cat $IN_DIR/$FILE | tr ' ' '\t' > $TMP.tab
cut -f 1 $TMP.tab > $TMP.1
cut -f 2 $TMP.tab > $TMP.2
$THIS/../trav $TMP.1 "echo -e '%f\t%h'" > $TMP.1h
$THIS/../trav $TMP.2 "echo -e '%f\t%h'" > $TMP.2h
paste $TMP.1h $TMP.2h > $TMP.ext
$THIS/../trav $TMP.ext "echo '$GENOME/%2/%1/%1.prot-univ $GENOME/%4/%3/%3.prot-univ'" > $TMP.filepair

$THIS/prots_pair2dissim  -log $LOG  -separate  $BLOSUM62  hmm-univ.list $TMP.filepair $TMP
cat $TMP | tr '\t' ' ' | sed 's/\.prot-univ / /g' | sed 's|'$GENOME'/[^/]\+/[^/]\+/||g' | tr ' ' '\t' > $OUT_DIR/$FILE


rm -f $LOG
rm $TMP*

