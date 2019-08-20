#!/bin/bash
if [ $# -ne 1 ]; then
  echo "Quality control of #1"
  echo "#1: incremental distance tree directory"
  echo "Requires: Time: O(n log^4(n))"
  exit 1
fi
INC=$1


