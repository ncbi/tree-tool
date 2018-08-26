#!/bin/csh -f

if ($# != 3) then
  echo "Print a .dm-file"
  echo "#1: File with numeric columns"
  echo "#2: Number precision"
  echo "#3: named(0/1)"
  exit 1
endif


set N = `wc -l $1`
if ($N[1] == 0) then
  echo "No data"
  exit 1
endif

set C = `head -1 $1 | wc -w`
@ C = $C[1] 
if ($3) then
  set NAME = name
  @ C = $C - 1
else
  set NAME = noname
endif

echo "ObjNum $N[1] $NAME nomult" 
echo "Attributes"
set i = 1
while ($i <= $C)
  echo "  V$i real $2"
  @ i = $i + 1
end
echo "Data"
cat $1

