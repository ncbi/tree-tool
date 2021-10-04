#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Convert NCBI taxonomy lineage to the 7-rank lineage"
  echo "#1: Input taxonomy feature file"
  exit 1
fi
IN=$1


grep '^01-00:\|^04-00:\|^07-00:\|^13-00:\|^18-00:\|^22-00:\|^26-00:' $IN
