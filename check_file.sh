#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: file path to check"
  exit 1
fi
F=$1


if [ ! -e $F ]; then  
  error "File or directory \"$F\" does not exist"
fi
