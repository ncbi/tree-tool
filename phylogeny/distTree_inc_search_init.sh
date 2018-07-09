#!/bin/csh -f

if ($# != 2) then
  echo "Initialize a new object for a distance tree"
  echo "Output (append): #1/hybrid"
  echo "#1: incremental distance tree directory"
  echo "#2: New object"
  exit 1
endif


set DIR = $1/search/$2


distTree_inc_new2hybrid.sh $1 $2 $DIR/hybrid $DIR/log
if ($?) exit 1
if (-e $DIR/log)  exit 1
if (-e $DIR/hybrid) then
  if (-z $DIR/hybrid) exit 1    
  flock $1/hybrid -c "cat $DIR/hybrid >> $1/hybrid"
  if ($?) exit 1
  rm -r $DIR/
  if ($?) exit 1
  exit 0
endif


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
