#!/bin/bash --noprofile
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "QC an object"
  echo "#1: object prefix"
  echo "#2: verbose: 0/1"
  exit 1
fi
PREF=$1
VERB=$2


if [ $VERB == 1 ]; then
  set -x
fi


NAME=$( basename $PREF )

ID==$( head -1 $PREF | sed 's/^>//1' | tr '\t' ' ' | cut -f 1 -d ' ' )
if [ "$NAME" != "$ID" ]; then
  error "File name '$NAME' must be the same as FASTA header identifier '$ID'"
fi
