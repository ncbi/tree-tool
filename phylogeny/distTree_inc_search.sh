#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Process a new object for a distance tree"
  echo "Update: append: #1/{leaf,dissim.add}"
  echo "        delete: #1/search/#2/"
  echo "        use: #1/log/#2"
  echo "#1: Directory containing search/"
  echo "#2: New object"
  echo "#3: Job number"
  echo "#4: Use grid (0/1)"
  exit 1
fi
INC=$1
OBJ=$2
JOB=$3
GRID=$4


#set -x


DIR=$INC/search/$OBJ


# Already done
if [ ! -e $DIR/request -a -s $DIR/dissim ]; then
  exit 0
fi

if [ ! -e $DIR/request ]; then
  error "$DIR/request does not exist" 
fi

if [ ! -e $DIR/dissim ]; then
  error "$DIR/dissim does not exist" 
fi


if [ -s $DIR/request ]; then
  cp /dev/null $INC/log/$OBJ
  if [ $GRID == 1 ]; then
    $QSUB_5 -N j$JOB "$THIS/distTree_inc_request.sh $INC $OBJ" > /dev/null  
      # Must be the last line of this script
  else
    $THIS/distTree_inc_request.sh $INC $OBJ
  fi
else
  # Finish
  cat $DIR/leaf >> $INC/leaf
  cat $DIR/dissim >> $INC/dissim.add
  rm -r $DIR/
fi

