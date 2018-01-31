#!/bin/csh -f

if ($# != 2) then
  echo "Build a distance tree incrementally"
  echo "Update: #1/"
  echo "#1: incremental distance tree directory"
  echo "#2: seed (>=1)"
  exit 1
endif



# Time: O(n log^3(n))
while (1)
  if (-e $1/stop) then
    echo "Stopped"
    exit 2
  endif
  
  set N = `ls $1/new/ | wc -l`
  echo "# To add: $N[1]  `date`  `date +%s`" >> $1/runlog  
  echo ""
  echo ""
  echo "# Total objects to add: $N[1] ..."
  if ($N[1] == 0) break  
  distTree_inc_new.sh $1 $2
  if ($?) exit 1
end
  

echo ""
echo "Complete optimization..."
# Time: O(n log^4(n))
makeDistTree  -data $1/  -optimize  -reinsert  -output_tree $1/tree.opt
if ($?) exit 1
