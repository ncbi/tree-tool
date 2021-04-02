#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Append inverted index of DNA mutations"
  echo "#1: input file with DNA mutations"
  echo "#2: directory with mutation files"
  exit 1
fi
IN=$1
DIR=$2


TMP=`mktemp`
#echo $TMP


NAME=`basename $IN`
grep -v '[^acgt0-9]' $IN > $TMP || true
$THIS/../trav $TMP -noprogress "echo $NAME >> $DIR/%f"


rm $TMP*



