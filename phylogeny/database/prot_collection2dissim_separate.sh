#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
GENOME="genome"
if [ $# -ne 5 ]; then
  echo 'Input: '$GENOME'/*/*/*.prot-univ, hmm-univ.list'
  echo "Output: #4/#1, #3.log/#1"
  echo "#1: file in #3 with pairs of genomes"
  echo "#2: 1 - BLOSUM62, 0 - PAM30"
  echo "#3: input directory"
  echo "#4: output directory"
  echo "#5: $GENOME/ is large (0/1)"
  exit 1
fi
FILE=$1
BLOSUM62_ARG=$2
IN_DIR=$3
OUT_DIR=$4
LARGE=$5


LOG=$IN_DIR.log/$FILE

BLOSUM62=""
if [ $BLOSUM62_ARG == 1 ]; then
  BLOSUM62="-blosum62"
fi


TMP=$( mktemp )
#echo $TMP
#set -x


# $TMP.filepair
if [ $LARGE == 1 ]; then
  < $IN_DIR/$FILE   tr ' ' '\t' > $TMP.tab
  cut -f 1 $TMP.tab > $TMP.1
  cut -f 2 $TMP.tab > $TMP.2
  $THIS/../../trav $TMP.1 "echo -e '%f\t%h'" -log $LOG > $TMP.1h
  $THIS/../../trav $TMP.2 "echo -e '%f\t%h'" -log $LOG > $TMP.2h
  paste $TMP.1h $TMP.2h > $TMP.ext
  $THIS/../../trav $TMP.ext "echo '$GENOME/%2/%1/%1.prot-univ $GENOME/%4/%3/%3.prot-univ'" -log $LOG > $TMP.filepair
else
  $THIS/../../trav $IN_DIR/$FILE "echo '$GENOME/%1/%1.prot-univ $GENOME/%2/%2.prot-univ'" -log $LOG > $TMP.filepair
fi

$THIS/../../dissim/prot_collection2dissim  -log $LOG  -separate  $BLOSUM62  hmm-univ.list $TMP.filepair $TMP

GENOME_SED='s|'$GENOME'/[^/]\+/||g'
if [ $LARGE == 1 ]; then
  GENOME_SED='s|'$GENOME'/[^/]\+/[^/]\+/||g'
fi
< $TMP  tr '\t' ' ' | sed 's/\.prot-univ / /g' | sed $GENOME_SED | tr ' ' '\t' > $OUT_DIR/$FILE


rm -f $LOG
rm $TMP*

