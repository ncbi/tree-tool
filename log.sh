#!/bin/bash
set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C

if [ $# -ne 1 ]; then
  echo "Create and populate #1.log/"
  echo "#1: Input directory with input files"
  echo "For a bash-only machine"
  exit 1
fi
DIR=$1


rm -rf $DIR.log
mkdir $DIR.log


TMP=`mktemp`

ls $DIR > $TMP
wc -l $TMP
while read ASM
do
  touch $DIR.log/$ASM
done < $TMP

rm -f $TMP*
