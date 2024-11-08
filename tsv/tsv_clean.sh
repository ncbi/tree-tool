#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Remove comment lines"
  echo "#1: tsv-file name"
  exit 1
fi
F=$1


TMP=$( mktemp )


grep '^#' $F | sort -u > $TMP
N=$( < $TMP wc -l )
if [ $N != 1 ]; then
  error "# Different headers: $N"
fi

cat $TMP
tail -n +2 $F | grep -v '^#' || true


rm $TMP
