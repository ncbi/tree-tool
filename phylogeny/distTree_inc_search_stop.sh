#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Stop processing a new object for a distance tree"
  echo "Update: append: #1/{leaf,dissim.add}"
  echo "        delete: #1/search/#2/"
  echo "#1: directory containing search/"
  echo "#2: new object"
  exit 1
fi
INC=$1
OBJ=$2


#set -x


DIR=$INC/search/$OBJ


if [ ! -e $DIR/dissim ]; then
  error "$DIR/dissim does not exist" 
fi

if [ ! -e $DIR/leaf ]; then
  error "$DIR/leaf does not exist" 
fi


cat $DIR/leaf >> $INC/leaf
cat $DIR/dissim >> $INC/dissim.add


rm -r $DIR/
