#!/bin/bash --noprofile
THIS=`dirname $0`
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "QC an object"
  echo "#1: file or directory with object data"
  exit 1
fi
FD=$1


NAME=`basename $FD`
PREF=$FD/$NAME

function check 
{
  SUF=$1
  #
  if [ ! -e $PREF.$SUF ]; then
    error "File $PREF.$SUF does not exist"
  fi
}

check "prot-univ"
check "prot"
check "univ"
