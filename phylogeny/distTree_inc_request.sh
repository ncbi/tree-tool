#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Process new objects for a distance tree"
  echo "Update: append: #1/search/#2/dissim"
  echo "                #1/log/#2"
  echo "        delete: #1/search/#2/request"
  echo "#1: incremental distance tree directory"
  echo "#2: new object"
  exit 1
fi
INC=$1
OBJ=$2


#set -x


DIR=$INC/search/$OBJ
LOG=$INC/log/$OBJ


$INC/pairs2dissim.sh $DIR/request "" $DIR/dissim.add $LOG &>> $LOG.pairs2dissim

# QC
N=`cat $DIR/request    | wc -l`
M=`cat $DIR/dissim.add | wc -l`
#if [ ! -s $DIR/dissim.add ]; then
if [ $N -ne $M ]; then
 #echo "Empty $DIR/dissim.add" > $LOG
  wc -l $DIR/request    >> $LOG
  wc -l $DIR/dissim.add >> $LOG
  exit 1
fi

rm $DIR/request

cat $DIR/dissim.add >> $DIR/dissim
rm $DIR/dissim.add


rm -f $LOG*
