#!/bin/csh -f

if ($# != 1) then
  echo "Process new objects for a distance tree: new/ -> leaf, dissim"
  echo "#1: distance tree data"
  exit 1
endif


set QC = ""  # -qc  
set RATE = 0.01   # PAR


set N = `ls $1/search/ | head -1`
if ("$N") then
  echo "$1/search/ is not empty"
  exit 1
endif

if (-e $1/dissim.add) then
  echo "$1/dissim.add exists"
  exit 1
endif
cp /dev/null $1/dissim.add
if ($?) exit 1

if (! -z $1/leaf) then
  echo "$1/leaf is not empty"
  exit 1
endif


echo "new/ -> search/ ..."

# Time: O(n) 
set OBJS = `grep -vc '^ *0x' $1/tree`
if ($?) exit 1
echo "# Objects: $OBJS"  

set INC = `echo "$OBJS * $RATE + 1" | bc -l | sed 's/\..*$//1'`  # PAR
echo "To add at this step: $INC"

ls $1/new/ > $1/new.list
if ($?) exit 1
set N = `wc -l $1/new.list`
if ($N[1] == 0) then
 #echo "New objects: $N[1]  Minimum: $INC"
  rm $1/new.list
  exit 
endif

setRandOrd $1/new.list 1 | head -$INC > $1/search.list
rm $1/new.list

trav -noprogress $1/search.list "mkdir $1/search/%f"
if ($?) exit 1
trav -noprogress $1/search.list "rm $1/new/%f"
if ($?) exit 1
rm $1/search.list


echo ""
echo "search/ -> leaf, dissim ..."


distTree_new $QC $1/ -init  
if ($?) exit 1
# Time of distTree_inc_request.sh: O(log(n)) per one new object
# Time: O(log^4(n)) per one new object, where n = # objects in the tree
set Iter = 0
while (1)
  set N = `ls $1/search/ | head -1`
  if ("$N" == "")  break  

  @ Iter = $Iter + 1
  echo ""
  echo "Iteration $Iter ..."
  
  set REQ_MAX = 2000  # PAR
  set REQ = `trav -noprogress $1/search "cat %d/%f/request" | wc -l`
  echo "# Requests: $REQ[1]"
  set GRID = 1
  if ($REQ[1] < $REQ_MAX)  set GRID = 0  

  rm -rf $1/log
  mkdir $1/log
  if ($?) exit 1

  trav $1/search "distTree_inc_search.sh $1 %f %n $GRID"
  if ($?) exit 1
  while ($GRID)
    sleep 10  # PAR
    set Q = `qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | head -1 | wc -l`
    if ($Q[1] == 0)  break
  end
  
  set L = `ls $1/log | wc -l`
  if ($L[1]) then
    echo "# Failed grid tasks: $L[1]"
    trav $1/log "distTree_inc_unsearch.sh $1 %f"
    if ($?) exit 1
  endif
  
  rm -r $1/log
  if ($?) exit 1
      
  distTree_new $QC $1/
  if ($?) exit 1
end


echo ""
echo "leaf, dissim.add -> tree, dissim ..."

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

#set VER = 161  

# Time: O(n log^2(n)) 
makeDistTree $QC  -data $1/  -remove_outliers $1/outlier.add  `cat $1/strong_outliers`  -reroot  -root_topological  -output_tree $1/tree.new  > $1/old/makeDistTree.$VER
if ($?) exit 1
mv $1/leaf $1/old/leaf.$VER
if ($?) exit 1
cp /dev/null $1/leaf
if ($?) exit 1
mv $1/tree.new $1/tree
if ($?) exit 1

trav -noprogress $1/outlier.add "cp /dev/null $1/outlier/%f"
if ($?) exit 1
rm $1/outlier.add
