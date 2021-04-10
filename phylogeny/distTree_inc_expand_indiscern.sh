#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Expand #1 by adding indiscernibile objects"
  echo "#1: incremental distance tree directory"
  echo "#2: sorted list of objects"
  exit 1
fi
INC=$1
LIST=$2


# QC
sort -c $LIST


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
cat $TMP.list


rm $TMP*
