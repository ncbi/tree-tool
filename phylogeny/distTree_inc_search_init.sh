#!/bin/csh -f

if ($# != 2) then
  echo "Initialize a new object for a distance tree"
  echo "#1: Directory containing search/"
  echo "#2: New object"
  exit 1
endif


set DIR = $1/search/$2


while (1) 
	$1/request_closest.sh $2 > $DIR/request
	if ($? == 0) break
	sleep 30
end

if (-z $DIR/request) then
  wc -l $DIR/request
  flock $1/alien -c "echo $2 >> $1/alien"
  if ($?) exit 1
  rm -r $DIR/
  if ($?) exit 1
else
	cp /dev/null $DIR/dissim
	if ($?) exit 1
endif
