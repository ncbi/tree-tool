#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "tar #1/ into #.tar.gz"
  echo "#1: directory"
  exit 1
fi
D=$1


tar -cvf $D.tar $D
gzip $D.tar

