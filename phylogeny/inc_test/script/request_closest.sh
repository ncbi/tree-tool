#!/bin/bash
THIS=`dirname $0`
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Find closest objects approximately"
  echo "#1: input object"
  exit 1
fi
OBJ=$1


# Fake procedure, only for a test
DT/phylogeny/tree2obj.sh $THIS/tree | sort -R | head -30 | sort | sed 's/^/'$OBJ' /1'

