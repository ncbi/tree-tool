#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print row fields in sepaeate rows with headers"
  echo "#1: tsv-file name"
  echo "#2: row key"
  exit 1
fi
F=$1
KEY="$2"


TMP=$( mktemp )


grep '^#' $F | tail -1 | sed 's/^#//1' | tr '\t' '\n' > $TMP.1
grep "$KEY" $F | tr '\t' '\n' > $TMP.2
N=$( < $TMP.2  wc -l )
if [$N -ne 1 ]; then
  cat $TMP.2
  error "Multiple rows"
fi
paste $TMP.1 $TMP.2


rm $TMP*