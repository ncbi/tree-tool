#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# != 5 ]; then
  echo "#1: Genome prefix"
  echo "#2: Genome-Hash database prefix"
  echo "#3: subset of Genome.id's (absolute pathname)"
  echo "#4: min. number of common hashes"
  echo "#5: output file (absolute pathname)"
  exit 1
fi
PREF=$1
DB=$2
SUBSET=$3
MIN_COMMON=$4
OUT=$5


TMP=$( mktemp )
#comment $TMP
#set -x


N=0
MAX=100   # PAR
i=0
touch $TMP.out
while [ $N -lt $MAX ]; do  
  i=$(( i + 1 ))
  case "$i" in
    "1") 
      SUF="CDS" 
      ;;
    "2") 
      SUF="PRT" 
      ;;
    "3") 
      SUF="HMM" 
      MIN_COMMON=0
      ;;
    *) 
      break 
      ;;
  esac
  if [ -e $DB.$SUF ]; then
    $THIS/../../objHash_find $DB.$SUF $PREF.hash-$SUF $MIN_COMMON $MAX  -superset $SUBSET >> $TMP.out
    sort -u $TMP.out > $TMP.out1
    mv $TMP.out1 $TMP.out
    N=$( < $TMP.out  wc -l )
  fi
done
cp $TMP.out $OUT


rm $TMP*
