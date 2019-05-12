#!/bin/bash
THIS=`dirname $0`
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: go"
  echo "Require: Time: O(n log^3(n))"
  exit 1
fi


TMP=`mktemp`
#echo $TMP


tree2obj.sh $THIS/tree > $TMP.list

grep "^>" $THIS/seq.fa | sed 's/^>//1' | sed 's/ .*$//1' | sort > $TMP.id

set +o errexit
uniq -c $TMP.id | grep -v '^ *1 ' > /dev/null
S=$?
set -o errexit
if [ $S -ne 1 ]; then
  exit 1
fi

diff $TMP.list $TMP.id


rm $TMP*


