#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
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
