#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: File to sort"
  exit 1
fi

TMP=`mktemp`
sort $1 > $TMP
mv $TMP $1
