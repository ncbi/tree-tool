#!/bin/bash
source /home/brovervv/code/cpp//bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: go"
  echo "Requires: Time: O(n log^3(n))"
  exit 1
fi


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
/home/brovervv/code/cpp/database/LocusQC.sh $INC $SERVER $DATABASE "Locus" "id" 4751 "28S"

