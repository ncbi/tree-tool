#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "#1: protein FASTA 1 (can be gzipped)"
  echo "#2: protein FASTA 2 (can be gzipped)"
  exit 1
fi
P1=$1
P2=$2


TMP=`mktemp`
#comment $TMP 


if file $P1 | grep -w "gzip" &> /dev/null ; then
  gunzip -c $P1 > $TMP.1
  P1=$TMP.1
fi

if file $P2 | grep -w "gzip" &> /dev/null ; then
  gunzip -c $P2 > $TMP.2
  P2=$TMP.2
fi

$THIS/../genetics/fasta2hash $P1 $TMP.h1  -prot_len_min 150  &> /dev/null 
$THIS/../genetics/fasta2hash $P2 $TMP.h2  -prot_len_min 150  &> /dev/null 

N1=`cat $TMP.h1 | wc -l`
N2=`cat $TMP.h2 | wc -l`
C=`$THIS/../setIntersect.sh $TMP.h1 $TMP.h2 1 | wc -l`

echo -e "$1\t$2\t$N1\t$N2\t$C"


rm $TMP*
