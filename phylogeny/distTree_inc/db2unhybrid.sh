#!/bin/bash
if [ $# -ne  1 ]; then
  echo "Print list of objects which have changed from hybrid to non-hybrid"
  echo "#1: output file"
  exit 1
fi
OUT=$1


exit 1