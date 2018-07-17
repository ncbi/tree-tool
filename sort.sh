#!/bin/csh -f

if ($# != 1) then
  echo "#1: File to sort"
  exit 1
endif

set TmpFNam = `mktemp`
sort $1 > $TmpFNam
if ($?) exit 1
mv $TmpFNam $1
if ($?) exit 1
