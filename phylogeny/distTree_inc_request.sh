#!/bin/bash --noprofile
THIS=$( dirname $0 )
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
INC=$( realpath $1 )
OBJ=$2


#set -x


DIR=$INC/search/$OBJ
LOG=$INC/log/$OBJ


if [ ! -e $LOG ]; then
  error "$LOG does not exist"
fi


$INC/pairs2dissim.sh $DIR/request "" $DIR/dissim.add $LOG &> $LOG.pairs2dissim


# QC
if [ ! -s $DIR/dissim.add ]; then
  error "Empty $DIR/dissim.add" >> $LOG
fi

N=$( cat $DIR/request    | wc -l )
M=$( cat $DIR/dissim.add | wc -l )
if [ $N -ne $M ]; then
  wc -l $DIR/request    >> $LOG
  wc -l $DIR/dissim.add >> $LOG
  exit 1
fi


rm $DIR/request

cat $DIR/dissim.add >> $DIR/dissim
rm $DIR/dissim.add


rm $LOG.pairs2dissim
rm -f $LOG
