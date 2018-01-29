#!/bin/csh -f

if ($# != 1) then
  echo "Print leaf statistics for an incremental distance tree"
  echo "#1: distance tree data"
  exit 1
endif

cat $1/hist/leaf* | sort -k 6 -g
if ($?) exit 1

