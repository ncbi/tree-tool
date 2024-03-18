#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
source $THIS/../qsub_env.sh
if [ $# -ne 4 ]; then
  echo "Add new objects in #1/new/ to a distance tree without the optimization of the original tree"
  echo "#1: incremental distance tree directory"
  echo "#2: input tree ('' = #1/tree)"
  echo "#3: output tree"
  echo "#4: output new leaf placement data file"
  echo "Time: O(n log^4(n))"
  exit 1
fi
INC=$1
IN="$2"
OUT=$3
PLACEMENT=$4


# Cf. distTree_inc_new.sh
GRID_MIN=`cat $INC/pairs2dissim.grid`
SEARCH_GRID_MIN=`cat $INC/object2closest.grid`  
QC=""  # -qc  
VARIANCE=`cat $INC/variance`

N=15
if [ -e $INC/threads ]; then
  N=`cat $INC/threads`
fi
THREADS="-threads $N"


N=`ls $INC/search/ | head -1`
if [ "$N" ]; then
  error "$INC/search/ is not empty"
fi

if [ -s $INC/leaf ]; then
  error "$INC/leaf is not empty"
fi

if [ -e $INC/dissim.add ]; then
  error "$INC/dissim.add exists"
fi


section "new/ -> search/"
$THIS/distTree_inc_new_list.sh $INC > $INC/search.list
OBJS=`cat $INC/search.list | wc -l`
wc -l $INC/search.list
if [ ! -s $INC/search.list ]; then
  warning "No objects to add"
  rm $INC/search.list
  exit 0
fi

$THIS/tree2obj.sh $INC/tree > $INC/tree.list

$THIS/../setIntersect.sh $INC/tree.list $INC/search.list 0 > $INC/tree-search.list
if [ -s $INC/tree-search.list ]; then
  wc -l $INC/tree-search.list
  error "Added objects are already in $INC/tree"
fi
rm $INC/tree-search.list


if [ $IN ]; then
  if [ -s $INC/tree.cur ]; then
    error "$INC/tree.cur exists"
  fi
  $THIS/tree2obj.sh $INC/tree > $INC/tree.list
  $THIS/tree2obj.sh $IN       > $INC/tree.list1
  if diff $INC/tree.list $INC/tree.list1; then
    rm $INC/tree.list
    rm $INC/tree.list1
    mv $INC/tree $INC/tree.cur
    cp $IN $INC/tree
  else
    error "$IN and $INC/tree have different objects"
  fi
fi


cp /dev/null $INC/dissim.add


$THIS/../trav  $THREADS  $INC/search.list "mkdir $INC/search/%f"
$THIS/distTree_inc_new_cmd.sh $INC "rm" $INC/search.list
rm $INC/search.list


section "search/ -> leaf, dissim"

N=`ls $INC/search/ | wc -l`
if [ $N -gt 0 ]; then
  rm -rf $INC/log/
  mkdir $INC/log
  $THIS/../trav  $INC/search "touch $INC/log/%f" 
  GRID=0
  if [ $N -lt $SEARCH_GRID_MIN ]; then
    $THIS/../trav  -step 1  $INC/search "$THIS/distTree_inc_search_init.sh $INC %f" 
  else
    GRID=1
    $THIS/../grid_wait.sh 1
    UL1=""
    if [ -e $INC/object2closest.sql ]; then
      UL1=",ul1=30"
    fi
    $THIS/../trav  -step 1  $INC/search "$QSUB_5$UL1  -N j%n  %Q$THIS/distTree_inc_search_init.sh $INC %f%Q > /dev/null" 
    WAIT=2000  # PAR
    if [ $SEARCH_GRID_MIN -le 200 ]; then  
      WAIT=7200
    fi
    $THIS/../qstat_wait.sh $WAIT 1
  fi
  $THIS/distTree_inc_new_log.sh $INC $GRID  
fi
rm $INC/tree.list


ITER=0
ITER_MAX=`echo $OBJS | awk '{printf "%d", 2*log($1)+3};'`
while [ $ITER -le $ITER_MAX ]; do
  # Time: O(log^4(n)) per one new object
  
  N=`ls $INC/search/ | wc -l`
  if [ $N == 0 ]; then
    break  
  fi

	ITER=$(( $ITER + 1 ))
  section "Iteration $ITER / $ITER_MAX"
  # use distTree_inc_request2dissim.sh ??
  REQ=`$THIS/../trav $INC/search "cat %d/%f/request" | wc -l`  
  echo "# Requests: $REQ"
  GRID=1
  if [ $REQ -lt $GRID_MIN ]; then
    GRID=0  
  fi

  rm -rf $INC/log/
  mkdir $INC/log

  if [ $GRID == 1 ]; then
	  $THIS/../grid_wait.sh 1
  fi
  $THIS/../trav  -step 1  $INC/search "$THIS/distTree_inc_search.sh $INC %f %n $GRID"
  WAIT=2000  # PAR
  if [ $GRID_MIN -lt 100 ]; then  
    WAIT=7200
  fi
  if [ $GRID == 1 ]; then
    $THIS/../qstat_wait.sh $WAIT 0
  fi
  
  $THIS/distTree_inc_new_log.sh $INC $GRID  

  echo "Processing new objects"
  $THIS/distTree_new $QC $INC/  -variance $VARIANCE
done
#$THIS/distTree_inc_dissim2indiscern.sh $INC $INC/dissim.add
rm $INC/dissim.add


section "leaf -> tree"
$THIS/makeDistTree $QC  $THREADS  -data $INC/  -variance $VARIANCE  -noqual  -output_tree $OUT
mv $INC/leaf $PLACEMENT
cp /dev/null $INC/leaf


if [ $IN ]; then
  mv $INC/tree.cur $INC/tree
fi


section "QC"
$INC/qc.sh 0
