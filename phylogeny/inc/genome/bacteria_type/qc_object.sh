#!/bin/bash --noprofile
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "QC an object"
  echo "#1: object prefix"
  exit 1
fi
PREF=$1


function check 
{
  local SUF=$1
  CPP_DIR/check_file.sh $PREF.$SUF 1
}

#check "hash-CDS"
#check "hash-HMM"
check "hash-PRT"
check "prot-univ"
check "prot.gz"
check "stat"
check "univ"
