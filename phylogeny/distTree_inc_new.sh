#!/bin/csh -f

if ($# != 2) then
  echo "Process new objects for a distance tree: new/ -> leaf, dissim"
  echo "exit 2: increments are completed"
  echo "#1: incremental distance tree directory"
  echo "#2: seed (>=1)"
  echo "Time: O(n log^4(n))"
  exit 1
endif


if ($2 <= 0)  exit 1


set GRID_MIN = `cat $1/grid_min`
set QC = ""  # -qc  
set RATE = 0.01   # PAR


if (1) then   
date
echo ""
top -b -n 1 | head -15
echo ""


set N = `ls $1/search/ | head -1`
if ("$N") then
  echo "$1/search/ is not empty"
  exit 1
endif

if (! -z $1/leaf) then
  echo "$1/leaf is not empty"
  exit 1
endif

if (-e $1/dissim.add) then
  echo "$1/dissim.add exists"
  exit 1
endif


set VER_OLD = `cat $1/version`
if ($?) exit 1

# Time: O(n) 
cp $1/tree $1/hist/tree.$VER_OLD
if ($?) exit 1

@ VER = $VER_OLD + 1
echo $VER > $1/version
if ($?) exit 1


echo "new/ -> search/ ..."

# Time: O(n) 
set OBJS = `grep -vc '^ *0x' $1/tree`
if ($?) exit 1
echo "# Objects: $OBJS"  

set INC = `echo "$OBJS * $RATE + 1" | bc -l | sed 's/\..*$//1'`  # PAR
echo "To add at this step: $INC"

ls $1/new/ > $1/new.list
if ($?) exit 1
if (-z $1/new.list) then
  rm $1/new.list
  exit 2
endif

cp /dev/null $1/dissim.add
if ($?) exit 1

setRandOrd $1/new.list  -seed $2 | head -$INC > $1/search.list
rm $1/new.list

trav -noprogress $1/search.list "mkdir $1/search/%f"
if ($?) exit 1
trav -noprogress $1/search.list "rm $1/new/%f"
if ($?) exit 1
rm $1/search.list


echo ""
echo "search/ -> leaf, dissim ..."

cp /dev/null $1/hybrid
if ($?) exit 1

set REQ = `ls $1/search | wc -l`
if ($REQ[1] > 20) then  # PAR
	trav  -step 1  $1/search "$QSUB_5,ul1=30  -N j%n  %QdistTree_inc_search_init.sh $1 %f%Q > /dev/null" 
	if ($?) exit 1  
	qstat_wait.sh 1
	if ($?) exit 1
else
	trav  -step 1  $1/search "distTree_inc_search_init.sh $1 %f"
	if ($?) exit 1  
endif

distTree_inc_hybrid.sh $1 $VER 
if ($?) exit 1
if (-e $1/hist/hybrid.$VER)  mv $1/hist/hybrid.$VER $1/hist/hybrid-new.$VER


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

wc -l $1/dissim.add
cat $1/dissim.add >> $1/dissim
if ($?) exit 1
rm $1/dissim.add
else
  set VER = 38  # PAR 
endif


# Time: O(n log^4(n)) 
set HYBRIDNESS_MIN = `cat $1/hybridness_min`
makeDistTree $QC  -threads 15 \
  -data $1/ \
  -optimize  -skip_len  -subgraph_iter_max 1 \
  -noqual \
  -delete_hybrids $1/hybrid  -delete_all_hybrids  -hybridness_min $HYBRIDNESS_MIN \
  -output_tree $1/tree.new \
  -dissim_request $1/dissim_request \
  > $1/hist/makeDistTree.$VER
if ($?) exit 1
  # -threads 20  # bad_alloc 
mv $1/leaf $1/hist/leaf.$VER
if ($?) exit 1
cp /dev/null $1/leaf
if ($?) exit 1
mv $1/tree.new $1/tree
if ($?) exit 1

cut -f 1 $1/hist/leaf.$VER | sort > $1/leaf.list
if ($?) exit 1
$1/objects_in_tree.sh $1/leaf.list 1
if ($?) exit 1
rm $1/leaf.list

echo ""
distTree_inc_hybrid.sh $1 $VER 
if ($?) exit 1

distTree_inc_unhybrid.sh $1 $VER 
if ($?) exit 1

echo ""
distTree_inc_request2dissim.sh $1 $1/dissim_request $1/dissim.add-req 
if ($?) exit 1
cat $1/dissim.add-req >> $1/dissim
if ($?) exit 1
rm $1/dissim.add-req
rm $1/dissim_request


if (-e $1/phen) then
  echo ""
  echo "Quality ..."
	tree2obj.sh $1/hist/tree.1 > $1/_init.list
	if ($?) exit 1	
	tree2obj.sh $1/tree > $1/_cur.list
	if ($?) exit 1	
	setMinus $1/_cur.list $1/_init.list >  $1/_delete.list
	if ($?) exit 1
	setMinus $1/_init.list $1/_cur.list >> $1/_delete.list
	if ($?) exit 1	
	rm $1/_init.list 
	rm $1/_cur.list
	makeDistTree  -input_tree $1/tree  -delete $1/_delete.list  -output_feature_tree $1/_feature_tree >& /dev/null
	if ($?) exit 1	
	rm $1/_delete.list
	makeFeatureTree  -input_tree $1/_feature_tree  -features $1/phen  -output_core $1/_core  -qual $1/_qual > $1/hist/makeFeatureTree.$VER
	if ($?) exit 1
	rm $1/_feature_tree
	rm $1/_core
	rm $1/_qual
	grep ' !' $1/hist/makeFeatureTree.$VER
	if ($?) exit 1
endif
