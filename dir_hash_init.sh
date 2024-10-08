#!/bin/bash --noprofile
THIS=`$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Create hashes in a directory for trav '%h'"
  echo "#1: directory"
  echo "#2: number of hashes (e.g., 1000)"
  exit 1
fi
DIR=$1
H=$2


$THIS/trav $H  -zero  -start 0  "mkdir $DIR/%n"  -threads 10  -step 1

