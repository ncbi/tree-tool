#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print total file size in M bytes"
  echo "#1: directory"
  exit 1
fi

du -s -m $1/*
