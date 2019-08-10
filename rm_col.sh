#!/bin/bash
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Remove columns in a tab-delimietd file"
  echo "#1: File"
  echo "#2: List of columns to be removed"
  exit 1
fi

TMP=`mktemp`
cut -f $2 --complement $1 > $TMP
mv $TMP $1
