#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Delete #1 without 'rm -r' which is slow"
  echo "#1: directory"
  exit 1
fi
DIR="$1"


if [ ! -d "$DIR" ]; then
  error "Directory '$DIR' does not exist"
fi


TMP=`mktemp -d`

# https://yonglhuang.com/rm-file/
rsync -a --delete $TMP/ "$DIR"
rmdir "$DIR"

rmdir $TMP


