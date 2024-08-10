#!/bin/bash --noprofile
THIS=`dirname $0`
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "QC an object"
  echo "#1: object path"
  exit 1
fi
F=$1


exit 1


function check 
{
  local SUF=$1
  if [ ! -e $F.$SUF ]; then
    error "File $F.$SUF does not exist"
  fi
}

#check "hash-CDS"
#check "hash-HMM"
#check "hash-PRT"
check "prot-univ"
check "prot.gz"
#check "stat"
check "univ"
