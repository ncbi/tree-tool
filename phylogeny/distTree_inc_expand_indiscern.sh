#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Expand #1 by adding indiscernibile objects"
  echo "#1: incremental distance tree directory" 
  echo "#2: sorted and distinct list of objects"
  echo "#3: restrict to the objects in #1/tree (0/1)"
  exit 1
fi
INC=$1
LIST=$2
IN_TREE=$3


# QC
sort -cu $LIST


TMP=`mktemp`
#echo $TMP
#set -x


awk '{print $1, $2};' $INC/indiscern >  $TMP
awk '{print $2, $1};' $INC/indiscern >> $TMP
sort -u $TMP | tr ' ' '\t' > $TMP.pair


cp $LIST $TMP.list
while true; do
  cp $TMP.list $TMP.expand
  join -1 1 -2 1 $TMP.list $TMP.pair | tr ' ' '\n' >> $TMP.expand
  sort -u $TMP.expand > $TMP.next
  if diff $TMP.list $TMP.next &> /dev/null ; then
    break
  fi
  mv $TMP.next $TMP.list
done

if [ $IN_TREE == 1 ]; then
  $THIS/tree2obj.sh $INC/tree > $TMP.tree
  $THIS/../setIntersect.sh $TMP.list $TMP.tree 0
else
  cat $TMP.list
fi


rm $TMP*
