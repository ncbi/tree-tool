#!/bin/csh -f

if ($# != 2) then
  echo "Restore #1/search/#2/"
  echo "#1: Directory containing search/"
  echo "#2: New object"
  exit 1
endif


set DIR = $1/search/$2


if (! -z $1/log/$2) then
  ls -laF $1/log/$2
  exit 1
endif

if (! -e $DIR/request) then
  echo "$DIR/request does not exist" 
  exit 1
endif

if (! -e $DIR/dissim) then
  echo "$DIR/dissim does not exist" 
  exit 1
endif


cp /dev/null $1/new/$2
if ($?) exit 1

rm -r $1/search/$2
if ($?) exit 1

