#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "List all directories in #1, recursively"
  echo "#1: directory"
  exit 1
fi
DIR=$1


echo $DIR
for D in `ls -dF $DIR/* | grep '/$' | sed 's|/$||1'`; do
  $THIS/ls_dirs.sh $D
done
