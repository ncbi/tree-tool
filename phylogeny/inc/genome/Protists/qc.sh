#!/bin/bash --noprofile
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: go"
  echo "Require: Time: O(n log^3(n))"
  exit 1
fi


INC=`dirname $0`


CPP_DIR/phylogeny/distTree_inc_indiscern_qc.sh $INC



exit 0

