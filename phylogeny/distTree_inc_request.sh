#!/bin/csh -f

if ($# != 2) then
  echo "Process new objects for a distance tree"
  echo "Update: append: #1/search/#2/dissim"
  echo "                #1/log/#2"
  echo "        delete: #1/search/#2/request"
  echo "#1: Directory containing search/"
  echo "#2: New object"
  exit 1
endif


set DIR = $1/search/$2
set LOG = $1/log/$2


$1/request2dissim.sh $DIR/request $DIR/dissim.add >> $LOG
if ($?) exit 1
if (-z $DIR/dissim.add) then
  echo "Empty $DIR/dissim.add" >> $LOG
  exit 1
endif

rm $DIR/request
if ($?) exit 1

cat $DIR/dissim.add >> $DIR/dissim
if ($?) exit 1
rm $DIR/dissim.add
if ($?) exit 1


rm $LOG
