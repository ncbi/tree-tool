#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Delete #1 without 'rm -r' which is slow"
  echo "#1: directory"
  exit 1
fi
DIR="$1"


if [ ! -e "$DIR" ]; then
  exit
fi


TMP=$( mktemp -d )

rsync -a --delete $TMP/ "$DIR"
rmdir "$DIR"

rmdir $TMP


