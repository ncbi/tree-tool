#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 1 ]; then
  echo "#1: qual-file of makeFeatureTree"
  exit 1
fi

grep ':' $1 | grep -v " +1 -0 /" | sed 's/:.*$//1' | sort | uniq -c
