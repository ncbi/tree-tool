#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "tsv_cat for a directory of .tsv-files"
  echo "#1: directory with .tsv-files"
  exit 1
fi
DIR=$1


TMP=`mktemp`


ls $DIR/* > $TMP.list
$THIS/tsv_cat $TMP.list "" -qc 


rm $TMP*

