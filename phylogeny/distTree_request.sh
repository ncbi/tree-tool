#!/bin/csh -f

if ($# != 2) then
  echo "Process new objects for a distance tree"
  echo "Update: append: #1/{leaf,dissim}"
  echo "        delete: #1/search/#2/"
  echo "#1: directory containing search/"
  echo "#2: new object"
  exit 1
endif


set DIR = $1/search/$2


if (! -e $DIR/request) then
  echo "$DIR/request does not exist"
  exit 1
endif

if (! -e $DIR/dissim) then
  echo "$DIR/dissim does not exist"
  exit 1
endif


if (-z $DIR/request) then
  cat $DIR/dissim >> $1/dissim.add
  if ($?) exit 1
  cat $DIR/leaf >> $1/leaf
  if ($?) exit 1
  rm -r $DIR/
  if ($?) exit 1
else
  mlst_request2dissim $DIR/request 28901 $DIR/dissim.add  # PAR
  if ($?) exit 1
  if (-z $DIR/dissim.add) then
    echo "Empty $DIR/dissim.add"
    exit 1
  endif

  rm $DIR/request
  if ($?) exit 1
  
  cat $DIR/dissim.add >> $DIR/dissim
  if ($?) exit 1
  rm $DIR/dissim.add
  if ($?) exit 1
endif



