#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Process new objects for a distance tree"
  echo "Update: append: #1/search/#2/dissim"
  echo "                #1/log/#2"
  echo "        delete: #1/search/#2/request"
  echo "#1: Directory containing search/"
  echo "#2: New object"
  exit 1
fi


#set -x


DIR=$1/search/$2
LOG=$1/log/$2


$1/request2dissim.sh $DIR/request $DIR/dissim.add $LOG &> $LOG.request2dissim
if [ ! -s $DIR/dissim.add ]; then
  echo "Empty $DIR/dissim.add" > $LOG
  exit 1
fi

rm $DIR/request

cat $DIR/dissim.add >> $DIR/dissim
rm $DIR/dissim.add


rm -f $LOG*
