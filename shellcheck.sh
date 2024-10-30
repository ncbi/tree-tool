#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Run shellcheck"
  echo "#1: bash script"
  exit 1
fi
S=$1


shellcheck  -e SC1008,SC1091,SC2010,SC2012,SC2045,SC2046,SC2086,SC2115,SC2129,SC2148,SC2155,SC2206,SC2207  -x  $S
  # SC2002,
