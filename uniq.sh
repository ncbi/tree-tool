#!/bin/bash
EXEC=~brovervv/code/cpp
source $EXEC/bash_common.sh
if [ $# != 1 ]; then
  echo "#1: File to sort and uniq"
  exit 1
fi

TMP=`mktemp`
cat $1 | tr '\t' ' ' | sed 's/ *$//1' | sed 's/  / /g' | sort | uniq | tr ' ' '\t' > $TMP
mv $TMP $1

