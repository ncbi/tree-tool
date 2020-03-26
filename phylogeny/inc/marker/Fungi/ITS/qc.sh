#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: go"
  echo "Requires: Time: O(n log^3(n))"
  exit 1
fi


INC=`dirname $0`
LocusQC.sh $INC Locus 4751 "5.8S"

