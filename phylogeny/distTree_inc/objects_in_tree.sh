#!/bin/bash --noprofile
#source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Record a list of objects in a database"
  echo "#1: list of objects"
  echo "#2: in_tree (0/null)"
  exit 1
fi
IN=$1
IN_TREE=$2


error "$0 is not implemented"
