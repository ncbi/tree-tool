#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print file names of empty files"
  echo "#1: directory"
  exit 1
fi
DIR=$1


#trav $DIR "if [ ! -s %d/%f ]; then echo %f; fi" -noprogress
ls -l $DIR/ | grep '  0 ' | awk '{print $9;}'

