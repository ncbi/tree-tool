#!/bin/tcsh -f

if ($# != 3) then
  echo "#1 - List1"
  echo "#2 - List2"
  echo "#3 - number(0/1)"
  echo "Make set-theoretic intersection of List1 and List2"
  exit 1
endif

set num = ""
if ($3) then
  set num = "-number"
endif


set TmpFNam=`mktemp` 

setMinus $num $1 $2 > $TmpFNam
if ($status) then
  cat $TmpFNam
  exit 1
endif

setMinus $num $1 $TmpFNam
if ($status) exit 1

rm $TmpFNam




