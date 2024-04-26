#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Difference of two tsv-files"
  echo "#1: tsv-file 1"
  echo "#2: tsv-file 2"
  echo "#3: key column names, must be common to #1 and #2 "
  exit 1
fi
F1=$1
F2=$2
KEYS="$3"


TMP=$( mktemp )
#comment $TMP
#set -x


$THIS/tsv2triple $F1  -keys "$KEYS"  -qc > $TMP.1
$THIS/tsv2triple $F2  -keys "$KEYS"  -qc > $TMP.2
diff $TMP.1 $TMP.2


rm $TMP*
