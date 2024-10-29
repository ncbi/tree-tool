#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "tar #1/ into #.tar.gz"
  echo "#1: directory"
  echo "#2: dereference (0/1)"
  exit 1
fi
D=$1
DEREF=$2


DEREF_PAR=""
if [ $DEREF == 1 ]; then
  DEREF_PAR="--dereference"
fi

tar -czf  $D.tar.gz $D  $DEREF_PAR


