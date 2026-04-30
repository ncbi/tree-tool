#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print <DNA1 name> <DNA2 name> <ANI>"
  echo "Not symmetric"
  echo "#1: DNA1 (FASTA)"
  echo "#2: DNA2 (FASTA)"
  exit 1
fi
DNA1=$1
DNA2=$2


TMP=$( mktemp )
#echo $TMP 


function len
{
  local F=$1
  $THIS/../genetics/fasta2len $F -noprogress | cut -f 2 | count | grep -w "^sum" | cut -f 2
}


ANI="nan"
$THIS/../genetics/dna_coverage.sh $DNA1 $DNA2 "combine" 0 100 '' 0 0 > $TMP
M=$( tail -n +2 $TMP | awk -F '\t' '{print $6 * $8};' | count | grep -w ^sum | cut -f 2 )
N=$( tail -n +2 $TMP | cut -f 6                       | count | grep -w ^sum | cut -f 2 )
if [ $N -gt 1000 ]; then  # PAR
  L1=$( len $DNA1 )
  L2=$( len $DNA2 )
  if [ $L1 -lt $L2 ]; then 
    L1=$L2
  fi
  GOOD=$( echo "$N / $L1 >= 0.7" | bc -l )  # PAR (Barbara)
  if [ $GOOD -eq 1 ]; then
    ANI=$( echo "scale=2; $M / $N" | bc )
  fi
fi

NAME1=$( basename $DNA1 )
NAME2=$( basename $DNA2 )
echo -e "$NAME1\t$NAME2\t$ANI"
 

rm $TMP*
