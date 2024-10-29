#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Remove a suffix from a file name"
  echo "#1: File name"
  echo "#2: File name suffix"
  exit 1
fi
FILE=$1
SUF=$2


DIR=$( dirname $FILE )
OUT=$( basename $FILE $SUF )
if [ "$FILE" != "$DIR/$OUT" ]; then
  mv $FILE $DIR/$OUT
fi

