#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print the list of objects in #1/new/"
  echo "#1: incremental distance tree directory"
  exit 1
fi
INC=$1


LARGE=`cat $INC/large`
if [ $LARGE == 1 ]; then
  $THIS/../trav $INC/new "ls %d/%f" | sort
else
  ls $INC/new/
fi
