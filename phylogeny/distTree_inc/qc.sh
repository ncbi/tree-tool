#!/bin/bash --noprofile
#source $BROVER_CPP/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality check"
  echo "#1: verbose (0/1)"
  echo "Requires: Time: O(n log^4(n))"
  exit 1
fi
VERB=$1


INC=$( dirname $0 )


TMP=$( mktemp )
if [ $VERB == 1 ]; then
  echo $TMP
  set -x
fi


if [ -e $INC/good ]; then
  sort -cu $INC/good
fi


tree2obj.sh $INC/tree > $TMP.tree
distTree_inc_new_list.sh $INC > $TMP.new
setIntersect.sh $TMP.tree $TMP.new 0 > $TMP.tree-new
if [ -s $TMP.tree-new ]; then
  wc -l $TMP.tree-new
  exit 1
fi


#warning "$0 is not implemented completely"


rm $TMP*
