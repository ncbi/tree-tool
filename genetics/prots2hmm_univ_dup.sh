#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find duplicate universal proteins"
  echo "#1: proteins"
  echo "#2: HMM library of universal proteins"
  echo "#3: number of cores for hmmsearch"
  echo "#4: output .tsv-file with dupliacate proteins"
  exit 1
fi
IN=$1
HMM_LIB=$2
CORES=$3
OUT=$4


TMP=$( mktemp )
#comment $TMP


hmmsearch  --tblout $TMP.hmmsearch  --domtblout $TMP.dom  --noali  -Z 10000  --cut_ga  --cpu $CORES  $HMM_LIB $IN &> $TMP.log
$THIS/hmmsearch2besthits $TMP.hmmsearch  -domtblout $TMP.dom  -duplicate $OUT > $TMP.annot


rm -fr $TMP*






