#!/bin/bash
set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C

if [ $# -ne 5 ]; then
  echo "Apply #3 to each file of #1, save in #2/, create #1.err"
  echo "#1: Directory with input files"
  echo "#2: Output directory"
  echo "#3: Log directory"
  echo "#4: 0 - #2/item is a file, 1 - #2/item is a directory"
  echo "#5: Executable with parameters: <Input> <Output> <Log file>"
  echo "For a bash-only machine"
  exit 1
fi
IN_DIR=$1
OUT_DIR=$2
LOG_DIR=$3
MAKE_DIR=$4
PROG=$5


TMP=`mktemp`


cp /dev/null $IN_DIR.err
ls $IN_DIR > $TMP
wc -l $TMP
while read ASM
do
  echo $ASM
  rm -rf $OUT_DIR/$ASM
  if [ $MAKE_DIR == 1 ]; then
    mkdir $OUT_DIR/$ASM
  fi
  set +o errexit
  $PROG $IN_DIR/$ASM $OUT_DIR/$ASM $LOG_DIR/$ASM
  S=$?
  set -o errexit
  if [ $S -ne 0 ]; then
    echo $ASM >> $IN_DIR.err
  fi
done < $TMP


rm -rf $TMP*

