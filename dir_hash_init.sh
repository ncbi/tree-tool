#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Create hashes in a directory for trav '%h'"
  echo "#1: directory"
  exit 1
fi
DIR=$1


$THIS/trav 1000 -zero -start 0 "mkdir $DIR/%n" -threads 10  -step 1

