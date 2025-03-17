#!/bin/bash --noprofile
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of #1"
  echo "#1: verbose (0/1)"
  echo "Requires: Time: O(n log^4(n))"
  exit 1
fi
VERB=$1


exit 0



