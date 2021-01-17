#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Execute a command for each object in a list of objects for #1/new/"
  echo "#1: incremental distance tree directory"
  echo "#2: command (e.g., touch, rm)"
  echo "#3: list of objects"
  exit 1
fi
INC=$1
CMD=$2
LIST=$3


H=""
LARGE=`cat $INC/large`
if [ $LARGE == 1 ]; then
  H="%h/"
fi

$THIS/../trav $LIST -threads 15 "$CMD $INC/new/$H%f"
