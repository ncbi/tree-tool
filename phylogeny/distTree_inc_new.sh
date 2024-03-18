#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
source $THIS/../qsub_env.sh
if [ $# -ne 1 ]; then
  echo "Add new objects to the distance tree: new/ -> leaf, dissim; optimize the tree"
  echo "#1: incremental distance tree directory"
  echo "Time: O(n log^4(n))"
  exit 1
fi
INC=$1



$THIS/../check_tmp.sh


# PAR
SEED=1  # >= 1


GRID_MIN=`cat $INC/pairs2dissim.grid`
SEARCH_GRID_MIN=`cat $INC/object2closest.grid`  
  # was: $GRID_MIN / 100
QC=""  # -qc  
RATE=0.015   # PAR
VARIANCE=`cat $INC/variance`
HYBRIDNESS_MIN=`cat $INC/hybridness_min`
REINSERT=""

$THIS/tree2obj.sh $INC/tree > $INC/tree.list

OBJS=`cat $INC/tree.list | wc -l`
ADD=`echo "$OBJS * $RATE" | bc -l | sed 's/\..*$//1'`  # PAR
if [ $ADD == 0 ]; then
  ADD=1
fi

N=15
if [ -e $INC/threads ]; then
  N=`cat $INC/threads`
fi
THREADS="-threads $N"


if true; then   
date
echo ""
echo ""


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


VER_OLD=`cat $INC/version`

# Time: O(n log(n)) 
cp $INC/tree $INC/hist/tree.$VER_OLD
if [ $VER_OLD != 1 ]; then
  gzip $INC/hist/tree.$VER_OLD
fi

VER=$(( $VER_OLD + 1 ))
echo $VER > $INC/version
super_section "version: $VER"


POS=$(( ${#VER} - 1 ))
if [ "${VER:$POS}" == "0" ]; then  # PAR
  REINSERT="-reinsert"
  warning "REINSERT\nTree quality temporarily decreases due to missing dissimilarities of close tree neighbors"
else
  section "new/ -> search/"
  echo "# Objects in tree: $OBJS"  
  echo "To add at this step: $ADD"

  $THIS/distTree_inc_new_list.sh $INC > $INC/new.list
  wc -l $INC/new.list

  cp /dev/null $INC/dissim.add

  $THIS/../setRandOrd $INC/new.list  -seed $SEED  -sigpipe | head -$ADD | sort > $INC/search.list
  wc -l $INC/search.list
  rm $INC/new.list

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
        # PAR
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


  ITER=0
  ITER_MAX=`echo $OBJS | awk '{printf "%d", log($1)+3};'`
  while [ $ITER -lt $ITER_MAX ]; do
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

    mkdir $INC/log

    if [ $GRID == 1 ]; then
  	  $THIS/../grid_wait.sh 1
    fi
    $THIS/../trav  -step 1  $INC/search "$THIS/distTree_inc_search.sh $INC %f %n $GRID"
    if [ $GRID == 1 ]; then
      WAIT=2000  # PAR
      if [ $GRID_MIN -le 200 ]; then  
        WAIT=7200
      fi
      $THIS/../qstat_wait.sh $WAIT 0
    fi
    
    $THIS/distTree_inc_new_log.sh $INC $GRID  

    $THIS/../trav  -step 1  $THREADS  $INC/search "$THIS/distTree_inc_search2bad.sh $INC %f"

    echo "Processing new objects"
    $THIS/distTree_new $QC $INC/  -variance $VARIANCE
  done
  $THIS/../trav  -step 1  $INC/search "$THIS/distTree_inc_search_stop.sh $INC %f"


  section "leaf, dissim.add -> tree, dissim"

  wc -l $INC/dissim.add
  cat $INC/dissim.add >> $INC/dissim
  $THIS/distTree_inc_dissim2indiscern.sh $INC $INC/dissim.add
  rm $INC/dissim.add
fi
else
  VER=`cat $INC/version`
fi
rm $INC/tree.list


DELETE_CRITERION_OUTLIERS=""
if [ -e $INC/delete_criterion_outliers ]; then
  DELETE_CRITERION_OUTLIERS="-delete_criterion_outliers $INC/outlier-criterion"
fi

HYBRID=""
DISSIM_BOUNDARY="NAN"
if [ "$HYBRIDNESS_MIN" != 0 ]; then
  DISSIM_BOUNDARY=`cat $INC/dissim_boundary`
	HYBRID="-hybrid_parent_pairs $INC/hybrid_parent_pairs  -delete_hybrids $INC/hybrid.new  -hybridness_min $HYBRIDNESS_MIN  -dissim_boundary $DISSIM_BOUNDARY"
  if [ ! -e $INC/delete_criterion_outliers ]; then
	  DELETE_CRITERION_OUTLIERS="-delete_criterion_outliers $INC/outlier-criterion  -criterion_outlier_num_max 1  -delete_deformation_outliers $INC/outlier-deformation  -deformation_outlier_num_max 1"
	fi
fi

DELETE=""
if [ -e $INC/outlier-genogroup ]; then
  wc -l $INC/outlier-genogroup
  DELETE="-delete $INC/outlier-genogroup  -check_delete"
fi

GOOD=""
if [ -e $INC/good ]; then
  $THIS/distTree_inc_expand_indiscern.sh $INC $INC/good 0 > $INC/good.expanded
  GOOD="-good $INC/good.expanded"
fi

# Time: O(n log^4(n)) 
$THIS/makeDistTree $QC  $THREADS  -data $INC/  -variance $VARIANCE \
  $DELETE \
  $REINSERT  -optimize  -skip_len  -subgraph_iter_max 2 \
  -noqual \
  $DELETE_CRITERION_OUTLIERS \
  $HYBRID \
  $GOOD \
  -output_tree     $INC/tree.new \
  -output_tree_tmp $INC/tree.tmp \
  -dissim_request $INC/dissim_request \
  > $INC/hist/makeDistTree.$VER
rm $INC/tree.tmp
mv $INC/leaf $INC/hist/leaf.$VER
cp /dev/null $INC/leaf
mv $INC/tree.new $INC/tree

if [ -s $INC/hist/leaf.$VER ]; then
  section "Database: new -> tree"
  cut -f 1 $INC/hist/leaf.$VER | sort > $INC/leaf.list
  $INC/objects_in_tree.sh $INC/leaf.list 1
  rm $INC/leaf.list
fi

if [ -e $INC/outlier-genogroup ]; then
  section "Database: genogroup outliers"
  wc -l $INC/outlier-genogroup
  $INC/objects_in_tree.sh $INC/outlier-genogroup null
  mv $INC/outlier-genogroup $INC/hist/outlier-genogroup.$VER
fi

if [ -e $INC/outlier-criterion ]; then
  section "Database: criterion outliers"
  wc -l $INC/outlier-criterion
  $INC/objects_in_tree.sh $INC/outlier-criterion null
  $THIS/../trav $INC/outlier-criterion "$INC/outlier2db.sh %f criterion"  
  mv $INC/outlier-criterion $INC/hist/outlier-criterion.$VER
fi

if [ -e $INC/outlier-deformation ]; then
  section "Database: deformation outliers"
  wc -l $INC/outlier-deformation
  cut -f 1 $INC/outlier-deformation > $INC/outlier-deformation.list
  $INC/objects_in_tree.sh $INC/outlier-deformation.list null
  $THIS/../trav $INC/outlier-deformation.list "$INC/outlier2db.sh %f deformation"  
  rm $INC/outlier-deformation.list
  mv $INC/outlier-deformation $INC/hist/outlier-deformation.$VER
fi

if [ "$HYBRIDNESS_MIN" != 0 ]; then
  section "Hybrid"
	$THIS/distTree_inc_hybrid.sh $INC 
 #echo "Unhybrid"
 #$THIS/distTree_inc_unhybrid.sh $INC 
fi

# Must be the last database change in this script
GENOGROUP_BARRIER=`cat $INC/genogroup_barrier`
if [ "$GENOGROUP_BARRIER" != "NAN" ]; then
  section "New genogroup outliers"
  if [ "$DISSIM_BOUNDARY" != "NAN" ]; then
    COMP=`echo "$DISSIM_BOUNDARY < $GENOGROUP_BARRIER" | bc`
    if [ $COMP == 1 ]; then
      error "dissim_boundary ($DISSIM_BOUNDARY) < genogroup_barrier ($GENOGROUP_BARRIER)"
    fi
  fi
  $THIS/tree2genogroup $INC/tree  $GENOGROUP_BARRIER  -genogroup_table $INC/genogroup_table
  $INC/genogroup2db.sh $INC/genogroup_table > $INC/outlier-genogroup  
  mv $INC/genogroup_table $INC/hist/genogroup_table.$VER
  gzip $INC/hist/genogroup_table.$VER
  if [ -s $INC/outlier-genogroup ]; then
    wc -l $INC/outlier-genogroup
  else
    rm $INC/outlier-genogroup
  fi
fi


section "QC $INC/indiscern"
$THIS/distTree_inc_indiscern_qc.sh $INC 0 0


section "Additional requests"
$THIS/distTree_inc_request2dissim.sh $INC $INC/dissim_request $INC/dissim.add-req
if [ -s $INC/dissim.add-req ]; then
  $THIS/distTree_inc_dissim2indiscern.sh $INC $INC/dissim.add-req
  cat $INC/dissim.add-req >> $INC/dissim
fi
rm $INC/dissim.add-req
rm $INC/dissim_request


$THIS/distTree_inc_tree1_quality.sh $INC


section "QC"
$INC/qc.sh 0


echo ""
NEW=`$THIS/distTree_inc_new_list.sh $INC | wc -l`
echo "# New objects: $NEW"
if [ $NEW == 0 -a -z "$REINSERT" ]; then
  touch $INC/finished
else
  rm -f $INC/finished
fi


