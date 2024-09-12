#!/bin/bash --noprofile
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "QC an object"
  echo "#1: object path"
  exit 1
fi
FD=$1


NAME=$( basename $FD )
PREF=$FD/$NAME

function check 
{
  local SUF=$1
  if [ ! -e $PREF.$SUF ]; then
    error "File $PREF.$SUF does not exist"
  fi
}

check "hash-CDS"
check "hash-HMM"
check "hash-PRT"
check "prot-univ"
check "prot.gz"
check "stat"
check "univ"
