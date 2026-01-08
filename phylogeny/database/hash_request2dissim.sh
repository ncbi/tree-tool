#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Add hash_request2dissim result to each line of #1"
  echo "#1: file with object pairs"
  echo "#2: directory with hashes <obj_hash>/<obj>/<obj>.<hash_type>"
  echo "#3: hash type" 
  echo "#4: intersection_min"
  echo "#5: ratio_min"
  exit 1
fi
PAIRS=$1
DIR=$2
TYPE=$3
INTERSECTION_MIN=$4
RATIO_MIN=$5


TMP=$( mktemp )
#comment $TMP


function pair2file
{
  local N=$1
  cut -f $N $PAIRS > $TMP.$N
  $THIS/../../file2hash $TMP.$N -file -append | awk '{printf "'$DIR'/%s/%s/%s.'$TYPE'\n", $1, $2, $2};' > $TMP.out.$N
}


pair2file 1
pair2file 2
paste $TMP.out.1 $TMP.out.2 > $TMP
hash_request2dissim $TMP $TMP.out  -intersection_min $INTERSECTION_MIN  -ratio_min $RATIO_MIN 
cut -f 3 $TMP.out > $TMP.dissim
paste $PAIRS $TMP.dissim


rm $TMP
