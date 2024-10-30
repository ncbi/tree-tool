#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Recompute hybridness"
  echo "#1: incremental distance tree directory"
  echo "#2: line of makeDistTree -delete_hybrids output file"
  exit 1
fi
INC=$1
LINE=( $2 )


TMP=$( mktemp )
#comment $TMP
#set -x


G1=${LINE[0]}
G2=${LINE[2]}
G3=${LINE[3]}

echo $G1 $G2 >  $TMP.req
echo $G2 $G3 >> $TMP.req
echo $G1 $G3 >> $TMP.req

$INC/pairs2dissim.sh $TMP.req "" $TMP.out $TMP.log &> /dev/null

L=( $( cut -f 3 $TMP.out | sort -n -r | tr '\n' ' ' ) )
R=$( echo "scale=2; ${L[0]} / (${L[1]} + ${L[2]})" | bc -l )
echo $G1 $G2 $G3 ${L[0]} ${L[1]} ${L[2]} $R | tr ' ' '\t'


rm $TMP*
