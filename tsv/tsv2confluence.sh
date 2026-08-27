#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "#1: tsv-file"
  echo "#2: replace underscores by spaces in the title (0/1)"
  exit 1
fi
F=$1
U2S=$2


TMP=$( mktemp )


sed 's/|/\\|/g' $F | tr '\t' '|' | sed 's/^/|/1' | sed 's/$/|/1' | sed 's/||/| |/g' | sed 's/||/| |/g' | sed 's/(\*)/(\\*)/g' > $TMP
head -1 $TMP | sed 's/^|#/|/1' | sed 's/|/||/g' > $TMP.head
if [ $U2S -eq 1 ]; then
  sed 's/_/ /g' $TMP.head
else
  cat $TMP.head
fi
tail -n +2 $TMP


rm $TMP