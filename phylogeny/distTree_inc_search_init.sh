#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Initialize a new object for a distance tree"
  echo "#1: incremental distance tree directory"
  echo "#2: New object"
  exit 1
fi


DIR=$1/search/$2


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
  flock $1/outlier-alien -c "echo $2 >> $1/outlier-alien"
  $1/outlier2db.sh $2 alien
  rm -r $DIR/
else
	cp /dev/null $DIR/dissim
fi
