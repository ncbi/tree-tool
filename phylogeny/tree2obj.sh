#!/bin/csh -f

if ($# != 1) then
  echo "#1: Distance tree directory"
  exit 1
endif

grep -v '^ *0x' $1/tree | sed 's/^ *//1' | sed 's/:.*$//1' | sort
if ($?) exit 1

