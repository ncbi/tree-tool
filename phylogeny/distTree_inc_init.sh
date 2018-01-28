#!/bin/csh -f

if ($# != 7) then
  echo "Initialize an incremental distance tree"
  echo "#1: Output directory"
  echo "#2: sorted list of objects"
  echo "#3: initial tree size"
  echo "#4: seed (>= 1)"
  echo "#5: request2dissim.sh line"
  echo "#6: objects_in_tree.sh line or ''"  
  echo "#7: request_closest.sh line or ''"
  echo "Time: 4 min. / 1000 objects"
  exit 1
endif


if (! -e $2 || -z $2) then
  echo "No list of objects"
  exit 1
endif
if ($3 <= 2)  exit 1
if ($4 <= 0)  exit 1
if ("$5" == "")  exit 1
if (("$6" == "") != ("$7" == ""))  exit 1


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

echo "$5" > $1/request2dissim.sh
if ($?) exit 1
chmod a+x $1/request2dissim.sh
if ($?) exit 1

if ("$6" != "") then
	echo "$6" > $1/objects_in_tree.sh
	if ($?) exit 1
	chmod a+x $1/objects_in_tree.sh
	if ($?) exit 1
endif

if ("$7" != "") then
	echo "$7" > $1/request_closest.sh
	if ($?) exit 1
	chmod a+x $1/request_closest.sh
	if ($?) exit 1
endif

echo "1" >> $1/version
if ($?) exit 1

mkdir $1/old
if ($?) exit 1

cp /dev/null $1/log
if ($?) exit 1


echo ""
echo "Initial tree ..."
setRandOrd $2 $4 | head -$3 | sort > $1/list.init

list2pairs $1/list.init > $1/dissim_request
if ($?) exit 1
wc -l $1/dissim_request

distTree_inc_request2dissim.sh $1 $1/dissim_request $1/dissim
if ($?) exit 1
rm $1/dissim_request

pairs2attr2 $1/dissim 1 cons 6 -distance > $1/data.dm
if ($?) exit 1

makeDistTree  -data $1/data  -dissim cons  -topology  -reroot  -root_topological  -remove_outliers $1/outlier.add  -output_tree $1/tree > $1/old/makeDistTree.1
if ($?) exit 1
trav -noprogress $1/outlier.add "cp /dev/null $1/outlier/%f"
if ($?) exit 1
rm $1/data.dm


if (-e $1/objects_in_tree.sh) then
	echo ""
	echo "Database ..."
	$1/objects_in_tree.sh $2 0
	if ($?) exit 1
	$1/objects_in_tree.sh $1/list.init 1
	if ($?) exit 1
	$1/objects_in_tree.sh $1/outlier.add 0
	if ($?) exit 1
	rm $1/outlier.add
endif


echo ""
echo "$1/new/ ..."
setMinus $2 $1/list.init > $1/new.list
if ($?) exit 1
rm $1/list.init
trav  -step 100  $1/new.list "cp /dev/null $1/new/%f"
if ($?) exit 1
rm $1/new.list


