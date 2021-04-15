#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: go"
  echo "Require: Time: O(n log^3(n))"
  exit 1
fi


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
CPP_DIR/phylogeny/database/LocusQC.sh $INC $SERVER $DATABASE Locus 2 "16S"
