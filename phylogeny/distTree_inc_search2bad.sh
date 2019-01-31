#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Remove a new contaminated object in a distance tree data structure indentified by inconsistent dissimilarities"
  echo "#1: incremental distance tree directory"
  echo "#2: New object"
  exit 1
fi


DIR=$1/search/$2


if [ ! -e $DIR/dissim ]; then
  exit 1
fi

awk '$3 == "nan"' $DIR/dissim > $DIR/dissim.bad

if [ -s $DIR/dissim.bad ]; then
  cat $DIR/dissim.bad >> $1/dissim.bad
  $1/outlier2db.sh $2 dissimilarity
  rm -r $DIR/
else
	rm $DIR/dissim.bad
fi
