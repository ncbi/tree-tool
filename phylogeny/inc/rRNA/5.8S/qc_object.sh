#!/bin/bash
THIS=`dirname $0`
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "QC an object"
  echo "#1: file or directory with object data"
  exit 1
fi
FD=$1


NAME=`basename $FD`

ID=`head -1 $FD | sed 's/^>//1' | tr '\t' ' ' | cut -f 1 -d ' '`
if [ "$NAME" != "$ID" ]; then
  error "File name '$NAME' must be the same as FASTA header identifier '$ID'"
fi
