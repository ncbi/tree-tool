#!/bin/bash
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print a row of a tab-delimnietd table where the top row is a header"
  echo "#1: table"
  echo "#2: row number (1-based)"
  exit 1
fi
TAB=$1
ROW=$2


L=`cat $TAB |wc -l`
L=$(( $L - 1 ))  # minus header
if [ $ROW -gt $L ]; then
  echo "Max. row = $L"
  exit 1
fi


TMP=`mktemp`
#echo $TMP


head -1 $TAB | sed 's/^#//1' | tr '\t' '\n' > $TMP.1

N=$(( $ROW + 1 ))
head -$N $TAB | tail -1 | tr '\t' '\n' > $TMP.2

paste $TMP.1 $TMP.2


rm $TMP*
