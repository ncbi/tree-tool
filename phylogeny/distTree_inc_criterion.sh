#!/bin/csh -f

if ($# != 1) then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  exit 1
endif

grep 'absCriterion =' -n $1/old/makeDistTree.* | sed 's|^'$1'/old/makeDistTree\.||1' | sed 's/:18:/ /1' | sort -n
if ($?) exit 1

