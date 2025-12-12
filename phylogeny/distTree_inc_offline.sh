#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Copy #1/ to #2/ and make #2/ empty and offline"
  echo "#1: input incremental distance tree directory"
  echo "#2: output incremental distance tree directory"
  exit 1
fi
FROM=$1
TO=$2


$THIS/../check_file.sh $FROM 0

if [ -e $TO ]; then
  error "TO exists"
fi
cp -r $FROM $TO

echo 1 > $TO/version

> $TO/tree
> $TO/dissim
> $TO/dissim.bad
> $TO/indiscern
> $TO/runlog
rm $TO/hist/*
rm -f $TO/tree.released

# Offline
> $TO/server
rm $TO/bulk
rm $TO/bulk_remote
rm $TO/database
rm $TO/objects_in_tree.sh
rm $TO/outlier2db.sh
rm $TO/genogroup2db.sh



