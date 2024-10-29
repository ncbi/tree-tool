#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Create #1.log/ and populate with items from #1/"
  echo "#1: directory with items"
  exit 1
fi
DIR=$1


mkdir $DIR.log
$THIS/trav $DIR  -step 10  -threads 10 "touch $DIR.log/%f"

