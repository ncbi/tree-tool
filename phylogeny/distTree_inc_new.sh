#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Process new objects for a distance tree: new/ -> leaf, dissim"
  echo "exit 2: increments are completed"
  echo "#1: incremental distance tree directory"
  echo "#2: seed (>=1)"
  echo "Time: O(n log^4(n))"
  exit 1
fi


if [ $2 -le 0 ]; then
  exit 1
fi


GRID_MIN=`cat $1/grid_min`
QC=""  # -qc  
RATE=0.01   # PAR


if [ 1 == 1 ]; then  # was: 1 == 1  
date
echo ""
top -b -n 1 | head -15
echo ""


N=`ls $1/search/ | head -1`
if [ $N ]; then
  echo "$1/search/ is not empty"
  exit 1
fi

if [ -s $1/leaf ]; then
  echo "$1/leaf is not empty"
  exit 1
fi

if [ -e $1/dissim.add ]; then
  echo "$1/dissim.add exists"
  exit 1
fi


VER_OLD=`cat $1/version`

# Time: O(n log(n)) 
cp $1/tree $1/hist/tree.$VER_OLD
if [ $VER_OLD != 1 ]; then
  gzip $1/hist/tree.$VER_OLD
fi

VER=$(( $VER_OLD + 1 ))
echo $VER > $1/version
echo "version: $VER"


echo "new/ -> search/ ..."

# Time: O(n) 
OBJS=`grep -vc '^ *0x' $1/tree`
echo "# Objects: $OBJS"  

INC=`echo "$OBJS * $RATE + 1" | bc -l | sed 's/\..*$//1'`  # PAR
echo "To add at this step: $INC"

ls $1/new/ > $1/new.list

cp /dev/null $1/dissim.add

setRandOrd $1/new.list  -seed $2  -sigpipe | head -$INC > $1/search.list
rm $1/new.list

trav -noprogress $1/search.list "mkdir $1/search/%f"
trav -noprogress $1/search.list "rm $1/new/%f"
rm $1/search.list


echo ""
echo "search/ -> leaf, dissim ..."

REQ=`ls $1/search | wc -l`
if [ $REQ -gt 20 ]; then  # PAR
	trav  -step 1  $1/search "$QSUB_5,ul1=30  -N j%n  %QdistTree_inc_search_init.sh $1 %f%Q > /dev/null" 
	qstat_wait.sh 1
else
	trav  -step 1  $1/search "distTree_inc_search_init.sh $1 %f"
fi

# Time: O(log^4(n)) per one new object, where n = # objects in the tree
ITER=0
while [ 1 == 1 ]; do
  N=`ls $1/search/ | wc -l`
  if [ $N == 0 ]; then
    break  
  fi

	ITER=$(( $ITER + 1 ))
  echo ""
  echo "Iteration $ITER ..."
  
  REQ=`trav -noprogress $1/search "cat %d/%f/request" | wc -l`
  echo "# Requests: $REQ"
  GRID=1
  if [ $REQ -lt $GRID_MIN ]; then
    GRID=0  
  fi

  rm -rf $1/log/
  mkdir $1/log

  trav  -step 1  $1/search "distTree_inc_search.sh $1 %f %n $GRID"
  if [ $GRID == 1 ]; then
    qstat_wait.sh 0
  fi
  
  L=`ls $1/log | wc -l`
  if [ $L -gt 0 ]; then
    echo "# Failed tasks: $L"
    if [ $GRID == 0 ]; then
      exit 1
    fi
	  # Try to fix grid problems
	  ls $1/log | grep -v '\.' > $1/log.list
    trav $1/log.list "distTree_inc_unsearch.sh $1 %f"
    rm $1/log.list
    trav $1/log "echo %d/%f; tail -20 %d/%f" > $1/log.out  # PAR
    head -21 $1/log.out # PAR
    rm $1/log.out
  fi
  
  rm -r $1/log/
      
  trav  -step 1  $1/search "distTree_inc_search2bad.sh $1 %f"

  echo "Processing new objects ..."
  distTree_new $QC $1/  -variance lin
done


echo ""
echo "leaf, dissim.add -> tree, dissim ..."

wc -l $1/dissim.add
cat $1/dissim.add >> $1/dissim
rm $1/dissim.add
else
  VER=`cat $1/version`
fi



DISSIM_BOUNDARY=`cat $1/dissim_boundary`

HYBRID=""
if [ -e $1/hybridness_min ]; then
	HYBRIDNESS_MIN=`cat $1/hybridness_min`
	HYBRID="-hybrid_parent_pairs $1/hybrid_parent_pairs  -delete_hybrids $1/hybrid.new  -delete_all_hybrids  -hybridness_min $HYBRIDNESS_MIN  -dissim_boundary $DISSIM_BOUNDARY"
fi

DELETE=""
if [ -e $1/outlier-genogroup ]; then
  wc -l $1/outlier-genogroup
  DELETE="-delete $1/outlier-genogroup  -check_delete"
fi

# Time: O(n log^4(n)) 
makeDistTree $QC  -threads 15  -data $1/  -variance lin \
  $DELETE \
  -optimize  -skip_len  -subgraph_iter_max 1 \
  -noqual \
  $HYBRID \
  -output_tree $1/tree.new \
  -dissim_request $1/dissim_request \
  > $1/hist/makeDistTree.$VER
  # -subgraph_iter_max 2
mv $1/leaf $1/hist/leaf.$VER
cp /dev/null $1/leaf
mv $1/tree.new $1/tree

echo ""
echo "Database: new -> tree ..."
cut -f 1 $1/hist/leaf.$VER | sort > $1/leaf.list
$1/objects_in_tree.sh $1/leaf.list 1
rm $1/leaf.list

if [ -e $1/outlier-genogroup ]; then
  echo "Database: genogroup outliers ..."
  $1/objects_in_tree.sh $1/outlier-genogroup null
  mv $1/outlier-genogroup $1/hist/outlier-genogroup.$VER
fi

echo ""
echo "Hybrid ..."
if [ -e $1/hybridness_min ]; then
	distTree_inc_hybrid.sh $1 
fi
echo "Unhybrid ..."
distTree_inc_unhybrid.sh $1 

# Must be the last database change in this script
GENOGROUP_BOUNDARY=`cat $1/genogroup_boundary`
if [ "$GENOGROUP_BOUNDARY" != "NAN" ]; then
  echo ""
  echo "New genogroup outliers ..."
  tree2species $1/tree  $GENOGROUP_BOUNDARY  -species_table $1/genogroup_table
  $1/genogroup2db.sh $1/genogroup_table > $1/outlier-genogroup  
  rm $1/genogroup_table 
  if [ -s $1/outlier-genogroup ]; then
    wc -l $1/outlier-genogroup
  else
    rm $1/outlier-genogroup
  fi
fi


echo ""
distTree_inc_request2dissim.sh $1 $1/dissim_request $1/dissim.add-req 
cat $1/dissim.add-req >> $1/dissim
rm $1/dissim.add-req
rm $1/dissim_request


distTree_inc_tree1_quality.sh $1


set +o errexit
L=`ls $1/new | head -1`
set -o errexit
if [ ! "$L" -a ! -e $1/outlier-genogroup ]; then
  exit 2
fi


