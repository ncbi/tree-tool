#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: Directory"
  exit 1
fi

du -s -m $1/*
