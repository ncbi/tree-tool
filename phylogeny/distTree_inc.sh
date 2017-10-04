#!/bin/csh -f

if ($# != 1) then
  echo "Build a distance tree incrementally"
  echo "Output: #1.log"
  echo "#1: distance tree data"
  exit 1
endif



rm -f $1.log
rm -f $1/leaf.old  # ??
while (1)
  sync
  set N = `ls $1/new/ | wc -l`
  echo "# New: $N[1]  `date`" >> $1.log  
  if ($N[1] == 0)  break
  
  echo ""
  echo ""
  echo ""
  echo "Adding $N[1] leaves ..."
  distTree_new.sh $1
  if ($?) exit 1
end
  
