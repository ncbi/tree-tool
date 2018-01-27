#!/bin/csh -f

if ($# != 2) then
  echo "Initialize a new object for a distance tree"
  echo "#1: Directory containing search/"
  echo "#2: New object"
  exit 1
endif


set DIR = $1/search/$2


if (-e $DIR/request) then
  echo "$DIR/request exists" 
  exit 1
endif

if (-e $DIR/dissim) then
  echo "$DIR/dissim exists" 
  exit 1
endif


$1/request_closest.sh $2 > $DIR/request
if ($?) exit 1
if (-z $DIR/request)  exit 1

cp /dev/null $DIR/dissim
if ($?) exit 1

