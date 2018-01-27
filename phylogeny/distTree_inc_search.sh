#!/bin/csh -f

if ($# != 4) then
  echo "Process a new object for a distance tree"
  echo "Update: append: #1/{leaf,dissim.add}"
  echo "        delete: #1/search/#2/"
  echo "#1: Directory containing search/"
  echo "#2: New object"
  echo "#3: Job number"
  echo "#4: Use grid (0/1)"
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
  cat $DIR/leaf >> $1/leaf
  if ($?) exit 1
  cat $DIR/dissim >> $1/dissim.add
  if ($?) exit 1
  rm -r $DIR/
  if ($?) exit 1
else
  cp /dev/null $1/log/$2
  if ($?) exit 1
  if ($4) then
    $QSUB -N j$3 "distTree_inc_request.sh $1 $2" > /dev/null  
  else
    distTree_inc_request.sh $1 $2
  endif
  if ($?) exit 1
endif

