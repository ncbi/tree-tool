#!/bin/csh -f

if ($# != 1) then
  echo "Process new objects for a distance tree: new/ -> leaf, dissim"
  echo "#1: distance tree data"
  echo "exit: 2 - too few new objects"
  exit 1
endif


set QC = ""  # -qc  
set RATE = 0.01   # PAR


echo "new/ -> search/ ..."

set N = `ls $1/search/ | head -1`
if ("$N") then
  echo "$1/search/ is not empty"
  exit 1
endif

# Time: O(n) 
set OBJS = `grep -vc '^ *0x' $1/tree`
if ($?) exit 1
echo "# Objects: $OBJS"  

set INC = `echo "$OBJS * $RATE + 1" | bc -l | sed 's/\..*$//1'`  # PAR
echo "To add: $INC"

ls $1/new/ > new.list
if ($?) exit 1
set N = `wc -l new.list`
if ($N[1] < $INC) then
  echo "New objects: $N[1]  Minimum: $INC"
  rm new.list
  exit 2
endif

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
# Time of distTree_inc_request.sh: O(log(n)) per one new object
# Time: O(log^4(n)) per one new object, where n = # leaves in the tree
set Iter = 0
while (1)
  set N = `ls $1/search/ | head -1`
  if ("$N" == "")  break  

  @ Iter = $Iter + 1
  echo ""
  echo "Iteration $Iter ..."
  
  mkdir $1/log
  if ($?) exit 1
  while (1)
    trav $1/search "distTree_inc_search.sh $1 %f %n"
    if ($?) exit 1
    while (1)
      sleep 10  # PAR
      set Q = `qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | head -1 | wc -l`
      if ($Q[1] == 0)  break
    end
    
    rmdir $1/log
    if ($? == 0) break
  end

  distTree_new $QC $1/
  if ($?) exit 1
end


echo ""
echo ""
echo "leaf, dissim -> tree, dissim ..."

set VER = `cat $1/version`
if ($?) exit 1

# Time: O(n log(n)) 
cp $1/dissim $1/dissim.old
if ($?) exit 1
# Time: O(n) 
cp $1/tree $1/old/tree.$VER
if ($?) exit 1

@ VER = $VER + 1
echo $VER > $1/version
if ($?) exit 1

wc -l $1/dissim.add
cat $1/dissim.add >> $1/dissim
if ($?) exit 1
rm $1/dissim.add

# Time: O(n log^2(n)) 
makeDistTree $QC  -data $1/  -remove_outliers $1/outlier.$VER  -output_tree $1/tree.new > $1/old/makeDistTree.$VER
if ($?) exit 1
mv $1/leaf $1/old/leaf.$VER
if ($?) exit 1
cp /dev/null $1/leaf
if ($?) exit 1
mv $1/tree.new $1/tree
if ($?) exit 1
cat $1/outlier.$VER >> $1/outlier
if ($?) exit 1
mv $1/outlier.$VER $1/old/
