#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "#1: File name"
  echo "#2: Line number"
  exit 1
fi
F=$1
L=$2


N=$( < $F wc -l )
if [ $L -gt $N ]; then
  error "Max. line = $N"
fi

head -$L $F | tail -1

