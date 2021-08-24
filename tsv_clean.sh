#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Remove comment lines"
  echo "#1: tsv-file name"
  exit 1
fi
F=$1


head -1 $F
tail -n +2 $F | grep -v '^#' || true
