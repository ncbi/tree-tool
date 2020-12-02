#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Initialize a new object for a distance tree"
  echo "#1: incremental distance tree directory"
  echo "#2: new object"
  exit 1
fi
INC=$1
OBJ=$2


DIR=$INC/search/$OBJ


while true; do
  set +o errexit
	$INC/request_closest.sh $OBJ | sed 's/$/\t'$OBJ'/1' > $DIR/request
  S=$?
  set -o errexit
	if [ $S == 0 ]; then
	  break
	fi
	sleep 30
done

if [ -s $DIR/request ]; then
	cp /dev/null $DIR/dissim
else
  wc -l $DIR/request
  $INC/outlier2db.sh $OBJ "alien"
  rm -r $DIR/
fi
