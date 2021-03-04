#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Symmetric best hits dissimilarity by k-mers: >= 0 or nan"
  echo "#1: genome/ (large directory)"
  echo "#2: Genome.id 1"
  echo "#3: Genome.id 2"
  exit 1
fi
GENOME=$1
ASM1=$2
ASM2=$3


TMP=`mktemp`
#echo $TMP 


function prepare
{
  ASM=$1
  OUT=$2
  #
  H=`$THIS/../file2hash $ASM`
  DIR=$GENOME/$H/$ASM
  gunzip $DIR/$ASM.prot.gz -c > $TMP.fa
  cut -f 1 $DIR/$ASM.univ > $TMP.univ
  $THIS/../genetics/extractFastaProt $TMP.fa $TMP.univ -remove > $OUT
}


prepare $ASM1 $TMP.1
prepare $ASM2 $TMP.2

$THIS/../dissim/symbet $TMP.1 $TMP.2  
  # -min_prot_len 150


rm $TMP*
