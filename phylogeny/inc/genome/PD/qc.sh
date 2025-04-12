#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
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


CPP_DIR/phylogeny/tree2obj.sh $INC/tree > $TMP.tree
CPP_DIR/phylogeny/distTree_inc_new_list.sh $INC > $TMP.new
CPP_DIR/setIntersect.sh $TMP.tree $TMP.new 0 > $TMP.tree-new
if [ -s $TMP.tree-new ]; then
  wc -l $TMP.tree-new
  exit 1
fi


#warning "$0 is not implemented completely"


rm $TMP*
