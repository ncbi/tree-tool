#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: Incremental distance tree directory"
  echo "Requires: Time: O(n log^4(n))"
  exit 1
fi
INC=$1


