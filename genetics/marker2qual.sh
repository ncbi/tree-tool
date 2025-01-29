#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print a good quality subset of eukaryotic marker proteins created by tblastn2marker_euk.sh"
  echo "#1: marker proteins (FASTA)"
  echo "#2: min. score to length ratio"
  echo "#3: output uniKernel file | ''"
  exit 1
fi
M=$1
T=$2
UNI=$3


TMP=$( mktemp )


grep '^>' $M | cut -f 1,7 -d ' '| sed 's/^>//1' | sed 's/ score=/\t/1' > $TMP.score
$THIS/fasta2len $M > $TMP.len
paste $TMP.len $TMP.score | awk -F '\t' '{OFS="\t"; print $1, $4/$2};' > $TMP.stat
if [ "$UNI" ]; then
  $THIS/../dm/conversion/cols2dm.sh $TMP.stat 0 5 1 > $TMP.dm
  $THIS/../dm/uniKernel $TMP "V2" -qc > $UNI
fi
awk '$2 > '$T $TMP.stat | cut -f 1 > $TMP.list
$THIS/filterFasta $M  -aa  -target $TMP.list  -len_min 20  -complexity_min 3


rm $TMP*  
