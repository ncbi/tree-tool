#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Restore #1/search/#2/"
  echo "#1: incremental distance tree directory"
  echo "#2: new object"
  exit 1
fi
INC=$1
OBJ=$2


H=""
if [ -e $INC/large ]; then
  H=`$THIS/../file2hash $OBJ`
  H="$H/"
fi
touch $INC/new/$H$OBJ

rm -r $INC/search/$OBJ/

