#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Initialize a new object for a distance tree"
  echo "#1: incremental distance tree directory"
  echo "#2: New object"
  exit 1
fi
INC=$1
OBJ=$2


DIR=$INC/search/$OBJ


while [ 1 == 1 ]; do
	set +o errexit
	$INC/request_closest.sh $OBJ > $DIR/request
  S=$?
	set -o errexit
	if [ $S == 0 ]; then
	  break
	fi
	sleep 30
done

if [ ! -s $DIR/request ]; then
  wc -l $DIR/request
  flock $INC/outlier-alien -c "echo $OBJ >> $INC/outlier-alien"
  $INC/outlier2db.sh $OBJ alien
  rm -r $DIR/
else
	cp /dev/null $DIR/dissim
fi
