#!/bin/csh -f

if ($# != 1) then
  echo "Build a distance tree incrementally"
  echo "Update: #1/"
  echo "Output: #1.log"
  echo "#1: distance tree data"
  exit 1
endif



rm -f $1.log
while (1)
  set N = `ls $1/new/ | wc -l`
  echo "# New left: $N[1]  `date`  `date +%s`" >> $1.log  
 #if ($N[1] == 0)  break
  
  echo ""
  echo ""
  echo ""
  echo "Adding total $N[1] leaves ..."
  distTree_inc_new.sh $1
  set S = $?
  if ($S == 2) break
  if ($S) exit 1
end
  
