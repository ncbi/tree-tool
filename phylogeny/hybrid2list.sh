#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Convert hybrid output of makeDistTree into a list and print"
  echo "#1: hybrid file"
  exit 1
fi
F=$1


TMP=$( mktemp )


awk '$7 == 1' $F | cut -f 1 >  $TMP
awk '$8 == 1' $F | cut -f 3 >> $TMP
awk '$9 == 1' $F | cut -f 4 >> $TMP
$THIS/../uniq.sh $TMP
cat $TMP


rm -rf $TMP*
