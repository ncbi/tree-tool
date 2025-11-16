#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: script"
  exit 1
fi
F=$1


grep '=$[1-9]' $F | sed 's/^.*=//1' | sort | uniq -c | grep -v '^ *1 ' || true
