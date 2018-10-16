#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Convert hybrid output of makeDistTree into a list and print"
  echo "#1: hybrid file"
  exit 1
fi


tmp=`mktemp`

cat $1 | awk '$7 == 1' | cut -f 1 >  $tmp
cat $1 | awk '$8 == 1' | cut -f 3 >> $tmp
cat $1 | awk '$9 == 1' | cut -f 4 >> $tmp
uniq.sh $tmp
cat $tmp

rm -rf $tmp*
