#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: tsv-file"
  exit 1
fi
F=$1


TMP=`mktemp`


cat $F | sed 's/|/\\|/g' | tr '\t' '|' | sed 's/^/|/1' | sed 's/$/|/1' | sed 's/||/| |/g' | sed 's/||/| |/g' > $TMP
head -1 $TMP | sed 's/|/||/g'
tail -n +2 $TMP


rm $TMP