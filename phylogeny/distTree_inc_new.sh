#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Process new objects for a distance tree: new/ -> leaf, dissim"
  echo "#1: Incremental distance tree directory"
  echo "Time: O(n log^4(n))+"
  exit 1
fi
INC=$1


# PAR
SEED=1  # >= 1


GRID_MIN=`cat $INC/grid_min`
QC=""  # -qc  
RATE=0.015   # PAR
VARIANCE=`cat $INC/variance`


if [ 1 == 1 ]; then  # was: 1 == 1  
date
echo ""
top -b -n 1 | head -15
echo ""


N=`ls $INC/search/ | head -1`
if [ $N ]; then
  echo "$INC/search/ is not empty"
  exit 1
fi

if [ -s $INC/leaf ]; then
  echo "$INC/leaf is not empty"
  exit 1
fi

if [ -e $INC/dissim.add ]; then
  echo "$INC/dissim.add exists"
  exit 1
fi


VER_OLD=`cat $INC/version`

# Time: O(n log(n)) 
cp $INC/tree $INC/hist/tree.$VER_OLD
if [ $VER_OLD != 1 ]; then
  gzip $INC/hist/tree.$VER_OLD
fi

VER=$(( $VER_OLD + 1 ))
echo $VER > $INC/version
echo "version: $VER"


echo "new/ -> search/ ..."

# Time: O(n) 
OBJS=`grep -vc '^ *0x' $INC/tree`
echo "# Objects: $OBJS"  

ADD=`echo "$OBJS * $RATE + 1" | bc -l | sed 's/\..*$//1'`  # PAR
echo "To add at this step: $ADD"

ls $INC/new/ > $INC/new.list

cp /dev/null $INC/dissim.add

$THIS/../setRandOrd $INC/new.list  -seed $SEED  -sigpipe | head -$ADD > $INC/search.list
rm $INC/new.list

$THIS/../trav -noprogress $INC/search.list "mkdir $INC/search/%f"
$THIS/../trav -noprogress $INC/search.list "rm $INC/new/%f"
rm $INC/search.list


echo ""
echo "search/ -> leaf, dissim ..."

# Time ??
REQ=`ls $INC/search | wc -l`
if [ $REQ -gt 20 ]; then  # PAR
	$THIS/../trav  -step 1  $INC/search "$QSUB_5,ul1=30  -N j%n  %Q$THIS/distTree_inc_search_init.sh $INC %f%Q > /dev/null" 
	$THIS/../qstat_wait.sh 2000 1
else
	$THIS/../trav  -step 1  $INC/search "$THIS/distTree_inc_search_init.sh $INC %f"
fi


# Time: O(log^3(n)) per one new object
# Limit by O(log(n)) iterations ??
ITER=0
while [ 1 == 1 ]; do
  N=`ls $INC/search/ | wc -l`
  if [ $N == 0 ]; then
    break  
  fi

	ITER=$(( $ITER + 1 ))
  echo ""
  echo "Iteration $ITER ..."
  
  REQ=`$THIS/../trav -noprogress $INC/search "cat %d/%f/request" | wc -l`
  echo "# Requests: $REQ"
  GRID=1
  if [ $REQ -lt $GRID_MIN ]; then
    GRID=0  
  fi

  rm -rf $INC/log/
  mkdir $INC/log

  $THIS/../trav  -step 1  $INC/search "$THIS/distTree_inc_search.sh $INC %f %n $GRID"
  if [ $GRID == 1 ]; then
    $THIS/../qstat_wait.sh 2000 0
  fi
  
  ls $INC/log | sed 's/\..*$//1' | sort | uniq > $INC/log.list
  L=`cat $INC/log.list | wc -l`
  if [ $L -gt 0 ]; then
    echo "# Failed tasks: $L"
    if [ $GRID == 0 ]; then
      exit 1
    fi
	  # Try to fix grid problems
    $THIS/../trav $INC/log.list "$THIS/distTree_inc_unsearch.sh $INC %f"
    $THIS/../trav $INC/log "echo ''; echo %d/%f; tail -20 %d/%f" > $INC/log.out  # PAR
    head -21 $INC/log.out # PAR
    rm $INC/log.out
  fi
  rm $INC/log.list
  
  rm -r $INC/log/
      
  $THIS/../trav  -step 1  $INC/search "$THIS/distTree_inc_search2bad.sh $INC %f"

  echo "Processing new objects ..."
  $THIS/distTree_new $QC $INC/  -variance $VARIANCE
done


echo ""
echo "leaf, dissim.add -> tree, dissim ..."

wc -l $INC/dissim.add
cat $INC/dissim.add >> $INC/dissim
rm $INC/dissim.add
else
  VER=`cat $INC/version`
fi


HYBRIDNESS_MIN=`cat $INC/hybridness_min`

HYBRID=""
if [ "$HYBRIDNESS_MIN" != 0 ]; then
  DISSIM_BOUNDARY=`cat $INC/dissim_boundary`
	HYBRID="-hybrid_parent_pairs $INC/hybrid_parent_pairs  -delete_hybrids $INC/hybrid.new  -hybridness_min $HYBRIDNESS_MIN  -dissim_boundary $DISSIM_BOUNDARY"
fi

DELETE=""
if [ -e $INC/outlier-genogroup ]; then
  wc -l $INC/outlier-genogroup
  DELETE="-delete $INC/outlier-genogroup  -check_delete"
fi

# Time: O(n log^3(n)) 
$THIS/makeDistTree $QC  -threads 15  -data $INC/  -variance $VARIANCE \
  $DELETE \
  -optimize  -skip_len  -subgraph_iter_max 2 \
  -noqual \
  $HYBRID \
  -delete_outliers $INC/outlier-criterion  -outlier_num_max 1 \
  -output_tree $INC/tree.new \
  -dissim_request $INC/dissim_request \
  > $INC/hist/makeDistTree.$VER
mv $INC/leaf $INC/hist/leaf.$VER
cp /dev/null $INC/leaf
mv $INC/tree.new $INC/tree

if [ -s $INC/hist/leaf.$VER ]; then
  echo ""
  echo "Database: new -> tree ..."
  cut -f 1 $INC/hist/leaf.$VER | sort > $INC/leaf.list
  $INC/objects_in_tree.sh $INC/leaf.list 1
  rm $INC/leaf.list
fi

if [ -e $INC/outlier-genogroup ]; then
  echo ""
  echo "Database: genogroup outliers ..."
  $INC/objects_in_tree.sh $INC/outlier-genogroup null
  mv $INC/outlier-genogroup $INC/hist/outlier-genogroup.$VER
fi

echo ""
echo "Database: criterion outlier ..."
wc -l $INC/outlier-criterion
$INC/objects_in_tree.sh $INC/outlier-criterion null
$THIS/../trav $INC/outlier-criterion "$INC/outlier2db.sh %f criterion"  
mv $INC/outlier-criterion $INC/hist/outlier-criterion.$VER

if [ "$HYBRIDNESS_MIN" != 0 ]; then
  echo ""
  echo "Hybrid ..."
	$THIS/distTree_inc_hybrid.sh $INC 
  echo "Unhybrid ..."
  $THIS/distTree_inc_unhybrid.sh $INC 
fi

# Must be the last database change in this script
GENOGROUP_BARRIER=`cat $INC/genogroup_barrier`
if [ "$GENOGROUP_BARRIER" != "NAN" ]; then
  echo ""
  echo "New genogroup outliers ..."
  $THIS/tree2genogroup $INC/tree  $GENOGROUP_BARRIER  -genogroup_table $INC/genogroup_table
  $INC/genogroup2db.sh $INC/genogroup_table > $INC/outlier-genogroup  
  mv $INC/genogroup_table $INC/hist/genogroup_table.$VER
  if [ -s $INC/outlier-genogroup ]; then
    wc -l $INC/outlier-genogroup
  else
    rm $INC/outlier-genogroup
  fi
fi


echo ""
$THIS/distTree_inc_request2dissim.sh $INC $INC/dissim_request $INC/dissim.add-req
if [ -s $INC/dissim.add-req ]; then
  grep -vwi nan $INC/dissim.add-req | grep -vwi inf >> $INC/dissim
fi
rm $INC/dissim.add-req
rm $INC/dissim_request


$THIS/distTree_inc_tree1_quality.sh $INC


NEW=`ls $INC/new | wc -l`
if [ $ADD -gt $NEW ]; then
  cp /dev/null $INC/finished
else
  rm -f $INC/finished
fi


