#!/bin/bash
set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C

if [ $# -ne 5 ]; then
  echo "Apply #3 to each file of #1, save in #2/"
  echo "#1: Directory with input files"
  echo "#2: Output directory"
  echo "#3: Log direcory"
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


rm -rf $OUT_DIR
mkdir $OUT_DIR


TMP=`mktemp`


ls $IN_DIR > $TMP
wc -l $TMP
N=0
M=0
while read ASM
do
  N=$(( $N + 1 ))
  M=$(( $M + 1 ))
  printf  "\r%d %s" $N $ASM
  rm -rf $OUT_DIR/$ASM
  if [ $MAKE_DIR == 1 ]; then
    mkdir $OUT_DIR/$ASM
  fi
  $PROG $IN_DIR/$ASM $OUT_DIR/$ASM $LOG_DIR/$ASM &
  if [ $M == 30 ]; then   # PAR
    wait
    M=0
  fi
done < $TMP
echo ""


rm -rf $TMP*

