#!/bin/bash --noprofile
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality check"
  echo "#1: go"
  echo "Requires: Time: O(n log^4(n))"
  exit 1
fi


warning "$0 is not implemented"
exit 0


INC=`dirname $0`


TMP=`mktemp`
#echo $TMP
#set -x


tree2obj.sh $INC/tree > $TMP.tree
ls $INC/new/ > $TMP.new
setIntersect.sh $TMP.tree $TMP.new > $TMP.tree-new
if [ -s $TMP.tree-new ]; then
  wc -l $TMP.tree-new
  exit 1
fi


rm $TMP*
