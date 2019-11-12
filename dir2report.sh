#!/bin/bash
set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C

if [ $# -ne 1 ]; then
  echo "Print a combined report"
  echo "#1: Input directory with reports"
  echo "Requries: each report file starts with a header line"
  echo "For a bash-only machine"
  exit 1
fi
DIR=$1


TMP=`mktemp`


ls $DIR > $TMP
if [ ! -s $TMP ]; then
  echo "No files"
  exit 1
fi
H=`head -1 $TMP`
head -1 $DIR/$H | sed -e 's/^/Item\t/1'
while read ASM
do
  tail -n +2 $DIR/$ASM | sed -e 's/^/'$ASM'\t/1'
done < $TMP


rm -rf $TMP*

