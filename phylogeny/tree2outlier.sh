#!/bin/csh -f

if ($# != 1) then
  echo "Print objects sorted by with criterion descending"
  echo "#1: Distance tree"
  exit 1
endif

grep 'leaf_error=' $1 | sed 's/:.* leaf_error=/#/1' | sed 's/^ *//1' | sed 's/ \+[^ ]*$//1' | tr '#' ' ' | sort -k 2 -r -g
if ($?) exit 1
