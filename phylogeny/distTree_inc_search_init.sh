#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Initialize a new object for a distance tree"
  echo "Output (append): #1/hybrid"
  echo "#1: incremental distance tree directory"
  echo "#2: New object"
  exit 1
fi


DIR=$1/search/$2


if [ 0 == 1 ]; then  # ??
	distTree_inc_new2hybrid.sh $1 $2 $DIR/hybrid $DIR/log 
	if [ -e $DIR/log ]; then
	  exit 1
	fi
	if [ -e $DIR/hybrid ]; then
	  if [ ! -s $DIR/hybrid ]; then
	    exit 1    
	  fi
	  flock $1/hybrid -c "cat $DIR/hybrid >> $1/hybrid"
	  rm -r $DIR/
	  exit 0
	fi
fi


while [ 1 == 1 ]; do
	set +o errexit
	$1/request_closest.sh $2 > $DIR/request
  S=$?
	set -o errexit
	if [ $S == 0 ]; then
	  break
	fi
	sleep 30
done

if [ ! -s $DIR/request ]; then
  wc -l $DIR/request
  flock $1/alien -c "echo $2 >> $1/alien"
  rm -r $DIR/
else
	cp /dev/null $DIR/dissim
fi
