#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Convert hybrid output of makeDistTree into a list and print"
  echo "#1: hybrid file"
  exit 1
fi


TMP=`mktemp`


cat $1 | awk '$7 == 1' | cut -f 1 >  $TMP
cat $1 | awk '$8 == 1' | cut -f 3 >> $TMP
cat $1 | awk '$9 == 1' | cut -f 4 >> $TMP
$THIS/../uniq.sh $TMP
cat $TMP


rm -rf $TMP*
