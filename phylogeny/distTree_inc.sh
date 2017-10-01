#!/bin/csh -f

if ($# != 1) then
  echo "Incrementally process new objects for a distance tree"
  echo "#1: distance tree data"
  exit 1
endif


set QC = -qc  # ""


distTree_new $QC $1/ -init
if ($?) exit 1
cp /dev/null $1/dissim.add
if ($?) exit 1
cp /dev/null $1/location.add
if ($?) exit 1


set Iter = 0
while (1)
  @ Iter = $Iter + 1
  echo ""
  echo ""
  echo "Iteration $Iter ..."
  
  set N = `ls $1/search | wc -l`
  if ($?) exit 1
  if ($N[1] == 0)  break  

  echo ""
  trav -step 1 $1/search "distTree_request.sh $1 %f"
  if ($?) exit 1
  echo ""

  distTree_new $QC $1/
  if ($?) exit 1
end
