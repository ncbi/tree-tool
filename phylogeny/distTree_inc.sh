#!/bin/csh -f

if ($# != 2) then
  echo "Build a distance tree incrementally"
  echo "Update: #1/"
  echo "#1: incremental distance tree directory"
  echo "#2: seed (>=1)"
  echo "Time: O(n log^5(n))"
  exit 1
endif



# Time: O(n log^4(n))
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
echo "Complete optimization ..."
echo "# Complete optimization  `date`  `date +%s`" >> $1/runlog  

set VER = `cat $1/version`
if ($?) exit 1

# Time: O(n) 
cp $1/tree $1/hist/tree.$VER
if ($?) exit 1

@ VER = $VER + 1
echo $VER > $1/version
if ($?) exit 1

set delete = ""
if (-e $1/delete) then
  wc -l $1/delete
  set delete = "-delete $1/delete"
endif

# Time: O(n log^5(n))
makeDistTree  -data $1/  $delete  -optimize  -reinsert  -output_tree $1/tree.new  -leaf_errors leaf_errors > $1/hist/makeDistTree.$VER
if ($?) exit 1
mv $1/tree.new $1/tree
if ($?) exit 1

if (-e $1/delete) then
  $1/objects_in_tree.sh $1/delete 0
	if ($?) exit 1
  mv $1/delete $1/hist/delete.$VER
	if ($?) exit 1
endif



if (-e $1/phen) then
	echo ""
	echo "Quality ..."
	tree_quality_phen.sh $1/tree $1/phen > $1/hist/tree_quality_phen.$VER
	if ($?) exit 1
	cat $1/hist/tree_quality_phen.$VER
	if ($?) exit 1
	set new_root = `grep '^New root: ' $1/hist/tree_quality_phen.$VER | sed 's/^New root: //1'`
	if ($?) exit 1
	set date = `date +%Y%m%d`
	if ($?) exit 1

	echo ""
	echo ""
	echo "New root: $new_root"
	echo ""
	echo ""
	makeDistTree  -data $1/  -reroot_at "$new_root"  -output_tree $1/tree.$date > /dev/null
	if ($?) exit 1
	echo ""
	tree_quality_phen.sh $1/tree.$date $1/phen > $1/hist/tree_quality_phen.$VER-rooted
	if ($?) exit 1
	cat $1/hist/tree_quality_phen.$VER-rooted
	if ($?) exit 1
	set new_root = `grep '^New root: ' $1/hist/tree_quality_phen.$VER-rooted | sed 's/^New root: //1'`
	if ("$new_root" != "") then
	  echo "Re-rooting must be idempotent"
	  exit 1
	endif
endif
