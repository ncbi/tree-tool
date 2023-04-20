#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Remove empty files"
  echo "#1: directory"
  exit 1
fi
DIR=$1


trav $DIR "if [ ! -s %d/%f ]; then rm %d/%f; fi" -noprogress
