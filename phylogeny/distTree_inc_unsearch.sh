#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Restore #1/search/#2/"
  echo "#1: Directory containing search/"
  echo "#2: New object"
  exit 1
fi


set DIR = $1/search/$2


if [ -s $1/log/$2 ]; then
  echo "Non-empty log:"
  ls -laF $1/log/$2
  exit 1
fi

if [ ! -e $DIR/dissim ]; then
  echo "$DIR/dissim does not exist" 
  exit 1
fi

cp /dev/null $1/new/$2
rm -r $1/search/$2

