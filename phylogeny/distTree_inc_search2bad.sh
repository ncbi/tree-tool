#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Remove a new contaminated object in a distance tree data structure indentified by inconsistent dissimilarities"
  echo "#1: incremental distance tree directory"
  echo "#2: new object"
  exit 1
fi
INC=$1
OBJ=$2


DIR=$INC/search/$OBJ


if [ ! -e $DIR/dissim ]; then
  error "$DIR/dissim does not exist"
fi

awk '$3 == "nan" || $3 == "-nan"' $DIR/dissim > $DIR/dissim.bad

if [ -s $DIR/dissim.bad ]; then
  flock $INC/dissim.bad -c "cat $DIR/dissim.bad >> $INC/dissim.bad"
  $INC/outlier2db.sh $OBJ 'incomparable'
  rm -r $DIR/
else
	rm $DIR/dissim.bad
fi
