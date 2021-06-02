#!/bin/bash --noprofile
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find closest objects"
  echo "#1: object"
  echo "#2: directory or ''"
  exit 1
fi
OBJ=$1
DIR="$2"


error "$0 is not implemented"

