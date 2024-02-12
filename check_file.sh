#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "#1: file or directory path to check"
  echo "#2: 1 - file, 0 - directory"
  exit 1
fi
F=$1
D=$2


if [  $D -eq 0 ]; then
  if [ ! -d $F ]; then  
    error "Directory \"$F\" does not exist"
  fi
else
  if [ ! -f $F ]; then  
    error "File \"$F\" does not exist"
  fi
fi
