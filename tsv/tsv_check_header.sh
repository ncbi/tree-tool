#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Check header"
  echo "#1: tsv-file name"
  echo "#2: Header line"
  exit 1
fi
F=$1
H="$2"


TMP=$( mktemp )


echo -e "$H" > $TMP.header
head -1 $F > $TMP.row1
diff $TMP.header $TMP.row1


rm $TMP*
