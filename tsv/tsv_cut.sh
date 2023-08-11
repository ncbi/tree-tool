#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Select rows by awk from a .tsv-file"
  echo "#1: .tsv-file name"
  echo "#2: comma-separated list of column numbers"
  echo "#3: make rows unique (0/1)"
  exit 1
fi
F=$1
C="$2"
UNIQ=$3


TMP=`mktemp`


AWKCOL=`echo $C | sed 's/,/,$/g' | sed 's/^/$/1'`

head -1 $F    | awk -F '\t' '{OFS="\t"; print '$AWKCOL'};' 
tail -n +2 $F | awk -F '\t' '{OFS="\t"; print '$AWKCOL'};' > $TMP
if [ $UNIQ == 1 ]; then
  sort -u $TMP
else
  cat $TMP
fi


rm $TMP*

