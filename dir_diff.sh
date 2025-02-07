#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
echo "Compare the files of two directories, except .-files"
  echo "#1: directory 1"
  echo "#2: directory 2"
  exit 1
fi
D1=$1
D2=$2


TMP=$( mktemp )


ls $D1 > $TMP.1
ls $D2 > $TMP.2

section "In $D1, but not in $D2"
$THIS/setMinus $TMP.1 $TMP.2 > $TMP.1-2
if [ -s $TMP.1-2 ]; then
  cat $TMP.1-2
fi

section "In $D2, but not in $D1"
$THIS/setMinus $TMP.2 $TMP.1 > $TMP.2-1
if [ -s $TMP.2-1 ]; then
  cat $TMP.2-1
fi


section "Comparing files"
$THIS/setIntersect.sh $TMP.1 $TMP.2 0 > $TMP.common
$THIS/trav -step 1  $TMP.common "diff $D1/%f $D2/%f || true"


rm $TMP*
