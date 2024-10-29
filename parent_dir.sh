#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print parent directory"
  exit 1
fi
F=$1


D=$( dirname $F )
P=$( realpath $D )
basename $P
