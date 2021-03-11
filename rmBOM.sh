#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Remove UTF-8 Byte Order Mark"
  echo "#1: Input text file"
  echo "#2: Output file"
  exit 1
fi
IN=$1
OUT=$2


if [ $IN == $OUT ];
  error "$IN = $OUT"
fi

cat $IN | sed 's/^\xEF\xBB\xBF//1' | sed 's/^\xFE\xFF//1' | sed 's/^\xFF\xFE//1' > $OUT
rm $IN
wc -l $OUT


