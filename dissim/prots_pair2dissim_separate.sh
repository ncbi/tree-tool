#!/bin/bash 
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo 'Input: genome/*/*.prot-univ, hmm-univ.list'
  echo "Output: #4/#1, #3.log/#1"
  echo "#1: File with pairs of genomes in #3"
  echo "#2: 1 - BLOSUM62, 0 - PAM30"
  echo "#3: Input directory"
  echo "#4: Output directory"
  exit 1
fi
FILE=$1
BLOSUM62_ARG=$2
IN_DIR=$3
OUT_DIR=$4


LOG=$IN_DIR.log/$FILE

BLOSUM62=""
if [ $BLOSUM62_ARG == 1 ]; then
  BLOSUM62="-blosum62"
fi


TMP=`mktemp`
#echo $TMP


rm -f $OUT_DIR/$FILE

cat $IN_DIR/$FILE | tr '\t' ' ' | awk -F ' ' '{printf "genome/%d/%d.prot-univ genome/%d/%d.prot-univ\n", $1, $1, $2, $2;}' > $TMP.filepair
$THIS/prots_pair2dissim  -log $LOG  -separate  $BLOSUM62  hmm-univ.list $TMP.filepair $TMP
cat $TMP | tr '\t' ' ' | sed 's/\.prot-univ / /g' | sed 's|genome/[^/]\+/||g' | tr ' ' '\t' > $OUT_DIR/$FILE


rm -f $LOG
rm $TMP*

