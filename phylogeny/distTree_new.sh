#!/bin/csh -f

if ($# != 1) then
  echo "Process new objects for a distance tree: new/ -> leaf, dissim"
  echo "#1: distance tree data"
  exit 1
endif


set QC = ""  # -qc  


echo "new/ -> search/ ..."

sync

set N = `ls $1/search/ | head -1`
if ("$N") then
  echo "$1/search/ is not empty"
  exit 1
endif

# Time: O(n) 
set OBJS = `grep -vc '^ *0x' subset500/tree`
if ($?) exit 1
echo "# Objects: $OBJS"  

set INC = `echo "$OBJS * 0.005 + 1" | bc -l | sed 's/\..*$//1'`  # PAR
echo "To add: $INC"

ls $1/new/ > new.list
if ($?) exit 1
setRandOrd new.list 1 | head -$INC > search.list
rm new.list

trav search.list "mkdir $1/search/%f"
if ($?) exit 1
trav search.list "rm $1/new/%f"
if ($?) exit 1
rm search.list


echo ""
echo ""
echo "search/ -> leaf, dissim ..."

if (-e $1/dissim.add) then
  echo "$1/dissim.add exists"
  exit 1
endif

if (! -z $1/leaf) then
  echo "$1/leaf is not empty"
  exit 1
endif

distTree_new $QC $1/ -init
if ($?) exit 1
# Time of distTree_request.sh: O(log(n)) per one new object
# Time: O(log^4(n)) per one new object
#   where n = # leaves in the tree
set Iter = 0
while (1)
  sync
  set N = `ls $1/search/ | head -1`
  if ("$N" == "")  break  

  @ Iter = $Iter + 1
  echo ""
  echo ""
  echo "Iteration $Iter ..."
  
  echo ""
  # qsub ??
  trav -step 1 $1/search "distTree_request.sh $1 %f"
  if ($?) exit 1
  
  echo ""
  distTree_new $QC $1/
  if ($?) exit 1
end


echo ""
echo ""
echo "leaf, dissim -> tree, dissim ..."

# Time: O(n log(n)) 
cp $1/dissim $1/dissim.old
if ($?) exit 1
# Time: O(n) 
cp $1/tree $1/tree.old
if ($?) exit 1

wc -l $1/dissim.add
cat $1/dissim.add >> $1/dissim
if ($?) exit 1
rm $1/dissim.add

# Time: O(n log^2(n)) 
makeDistTree $QC  -data $1/  -output_tree $1/tree.new
if ($?) exit 1
echo ""     >> $1/leaf.old  # ??
cat $1/leaf >> $1/leaf.old  # ??
if ($?) exit 1
cp /dev/null $1/leaf
if ($?) exit 1
mv $1/tree.new $1/tree
if ($?) exit 1
