#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print the list of objects in #1/new/"
  echo "#1: incremental distance tree directory"
  exit 1
fi
INC=$1


if [ -e $INC/large ]; then
  N=$( ls $INC/new | wc -l )
  if [ $N -ne 1000 ]; then
    error "$INC/new/ has no 1000 subdirectories (it has $N)"
  fi
  $THIS/../trav $INC/new "ls %d/%f" | sort
else
  ls $INC/new/
fi
