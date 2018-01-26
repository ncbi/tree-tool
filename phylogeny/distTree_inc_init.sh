#!/bin/csh -f

if ($# != 5) then
  echo "Initialize an incremental distance tree"
  echo "#1: Output distance tree data directory"
  echo "#2: request2dissim.sh line"
  echo "#3: sorted list of objects"
  echo "#4: initial tree size"
  echo "#5: seed (>= 1)"
  exit 1
endif


if (! -e $3 || -z $3) then
  echo "No list of objects"
  exit 1
endif
if ($4 <= 2)  exit 1
if ($5 <= 0)  exit 1


mkdir $1
if ($?) exit 1


cp /dev/null $1/tree
if ($?) exit 1

mkdir $1/outlier
if ($?) exit 1

cp /dev/null $1/dissim
if ($?) exit 1

cp /dev/null $1/leaf
if ($?) exit 1

mkdir $1/new
if ($?) exit 1

mkdir $1/search
if ($?) exit 1

echo "$2" > $1/request2dissim.sh
if ($?) exit 1
chmod a+x $1/request2dissim.sh
if ($?) exit 1

echo "1" >> $1/version
if ($?) exit 1

mkdir $1/old
if ($?) exit 1


# initial tree 
setRandOrd $3 $5 | head -$4 | sort > $1/list.init

list2pairs $1/list.init > $1/dissim_request
if ($?) exit 1
wc -l $1/dissim_request

distTree_inc_dissim_request.sh $1 $1/dissim_request $1/dissim
if ($?) exit 1
rm $1/dissim_request

pairs2attr2 $1/dissim 1 cons 6 -distance > $1/data.dm
if ($?) exit 1

makeDistTree  -data $1/data  -dissim cons  -topology  -reroot  -root_topological  -remove_outliers $1/outlier.add  -output_tree $1/tree > $1/old/makeDistTree.1
if ($?) exit 1
trav -noprogress $1/outlier.add "cp /dev/null $1/outlier/%f"
if ($?) exit 1
rm $1/outlier.add
rm $1/data.dm


setMinus $3 $1/list.init > $1/new.list
if ($?) exit 1
rm $1/list.init
trav  -step 100  $1/new.list "cp /dev/null $1/new/%f"
if ($?) exit 1
rm $1/new.list


