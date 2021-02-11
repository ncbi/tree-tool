#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Symmetric best hits dissimilarity by 5-mers:>= 0 or nan"
  echo "#1: gzip'ed protein FASTA 1"
  echo "#2: gzip'ed protein FASTA 2"
  exit 1
fi
PROTGZ1=$1
PROTGZ2=$2


TMP=`mktemp`
#echo $TMP 


gunzip $PROTGZ1 -c > $TMP.1
gunzip $PROTGZ2 -c > $TMP.2

$THIS/symbet $TMP.1 $TMP.2 


rm $TMP*
