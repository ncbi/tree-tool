#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "diff on selected columns of two .tsv-files"
  echo "#1: .tsv-file name #1"
  echo "#2: .tsv-file name #2"
  echo "#3: cut parameters"
  exit 1
fi
F1=$1
F2=$2
CUT="$3"


TMP=$( mktemp )
#comment $TMP


cut -f $CUT $F1 > $TMP.1
cut -f $CUT $F2 > $TMP.2
differ $TMP.1 $TMP.2


rm $TMP*

