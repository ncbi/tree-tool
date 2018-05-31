#!/bin/csh -f

if ($# != 2) then
  echo "Process new objects for a distance tree: new/ -> leaf, dissim"
  echo "#1: incremental distance tree directory"
  echo "#2: seed (>=1)"
  echo "Time: O(n log^2(n))"
  echo "Cumulative time: O(n log^4(n))"
  exit 1
endif


if ($2 <= 0)  exit 1


set GRID_MIN = `cat $1/grid_min`
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

setRandOrd $1/new.list $2 | head -$INC > $1/search.list
rm $1/new.list

trav -noprogress $1/search.list "mkdir $1/search/%f"
if ($?) exit 1
trav -noprogress $1/search.list "rm $1/new/%f"
if ($?) exit 1
rm $1/search.list


echo ""
echo "search/ -> leaf, dissim ..."

set REQ = `ls $1/search | wc -l`
if ($REQ[1] > 200) then  # PAR
	trav  -step 1  $1/search "$QSUB_5 -N j%n %QdistTree_inc_search_init.sh $1 %f%Q > /dev/null" 
	if ($?) exit 1  
	qstat_wait.sh 1
	if ($?) exit 1
else
	trav  -step 1  $1/search "distTree_inc_search_init.sh $1 %f"
	if ($?) exit 1  
endif


# Time: O(log^4(n)) per one new object, where n = # objects in the tree
set Iter = 0
while (1)
  set N = `ls $1/search/ | head -1`
  if ("$N" == "")  break  

  @ Iter = $Iter + 1
  echo ""
  echo "Iteration $Iter ..."
  
  set REQ = `trav -noprogress $1/search "cat %d/%f/request" | wc -l`
  echo "# Requests: $REQ[1]"
  set GRID = 1
  if ($REQ[1] < $GRID_MIN)  set GRID = 0  

  rm -rf $1/log/
  mkdir $1/log
  if ($?) exit 1

  trav  -step 1  $1/search "distTree_inc_search.sh $1 %f %n $GRID"
  if ($?) exit 1
  if ($GRID) then
    qstat_wait.sh 0
	  if ($?) exit 1
  endif
  
  set L = `ls $1/log | wc -l`
  if ($L[1]) then
    echo "# Failed tasks: $L[1]"
   #exit 2
    if (! $GRID) exit 1
	  # Try to fix grid problems
    trav $1/log "distTree_inc_unsearch.sh $1 %f"
    if ($?) exit 1
    trav $1/log "echo %d/%f; tail -20 %d/%f" | head -21
  endif
  
  rm -r $1/log/
  if ($?) exit 1
      
  echo "Processing new objects ..."
  distTree_new $QC $1/
  if ($?) exit 1
end


echo ""
echo "leaf, dissim.add -> tree, dissim ..."

set VER = `cat $1/version`
if ($?) exit 1

# Time: O(n) 
cp $1/tree $1/hist/tree.$VER
if ($?) exit 1

@ VER = $VER + 1
echo $VER > $1/version
if ($?) exit 1

wc -l $1/dissim.add
cat $1/dissim.add >> $1/dissim
if ($?) exit 1
rm $1/dissim.add

#set VER = 254


# Time: O(n log^2(n)) 
# -profile
makeDistTree $QC  -data $1/ \
  -optimize  -skip_len  -subgraph_fast  -max_subgraph_iter 2 \
  -reroot  -root_topological \
  -noqual \
  -output_tree $1/tree.new \
  -dissim_request $1/dissim_request \
  > $1/hist/makeDistTree.$VER
if ($?) exit 1
  # -remove_outliers $1/outlier.add \
mv $1/leaf $1/hist/leaf.$VER
if ($?) exit 1
cp /dev/null $1/leaf
if ($?) exit 1
mv $1/tree.new $1/tree
if ($?) exit 1

echo ""
cut -f 1 $1/hist/leaf.$VER > $1/leaf.list
if ($?) exit 1
$1/objects_in_tree.sh $1/leaf.list 1
if ($?) exit 1
rm $1/leaf.list
#$1/objects_in_tree.sh $1/outlier.add 0
#if ($?) exit 1

if (0) then
	echo ""
	trav  -noprogress  $1/outlier.add "cp /dev/null $1/outlier/%f"
	if ($?) exit 1
	rm $1/outlier.add
endif

distTree_inc_request2dissim.sh $1 $1/dissim_request $1/dissim.add-req 
if ($?) exit 1
cat $1/dissim.add-req >> $1/dissim
if ($?) exit 1
rm $1/dissim.add-req
rm $1/dissim_request


