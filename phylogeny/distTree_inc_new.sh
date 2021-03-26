#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Process new objects for a distance tree: new/ -> leaf, dissim"
  echo "#1: incremental distance tree directory"
  echo "Time: O(n log^4(n))"
  exit 1
fi
INC=$1



$THIS/../check_tmp.sh


# PAR
SEED=1  # >= 1


GRID_MIN=`cat $INC/grid_min`
QC=""  # -qc  
RATE=0.015   # PAR
VARIANCE=`cat $INC/variance`
OBJS=`$THIS/tree2obj.sh $INC/tree | wc -l`
ADD=`echo "$OBJS * $RATE" | bc -l | sed 's/\..*$//1'`  # PAR
if [ $ADD == 0 ]; then
  ADD=1
fi


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


section "new/ -> search/"
echo "# Objects in tree: $OBJS"  
echo "To add at this step: $ADD"

$THIS/distTree_inc_new_list.sh $INC > $INC/new.list
wc -l $INC/new.list

cp /dev/null $INC/dissim.add

$THIS/../setRandOrd $INC/new.list  -seed $SEED  -sigpipe | head -$ADD | sort > $INC/search.list
wc -l $INC/search.list
rm $INC/new.list

$THIS/../trav  -threads 15  $INC/search.list "mkdir $INC/search/%f"
$THIS/distTree_inc_new_cmd.sh $INC "rm" $INC/search.list
rm $INC/search.list


section "search/ -> leaf, dissim"

N=`ls $INC/search/ | wc -l`
if [ $N -gt 0 ]; then
  rm -rf $INC/log/
  mkdir $INC/log
  $THIS/../trav  $INC/search "touch $INC/log/%f" 
  SEARCH_GRID_MIN=$(( $GRID_MIN / 100 ))  # PAR
  GRID=0
  if [ $N -lt $SEARCH_GRID_MIN ]; then
    $THIS/../trav  -step 1  $INC/search "$THIS/distTree_inc_search_init.sh $INC %f" 
  else
    GRID=1
    $THIS/../grid_wait.sh 1
    UL1=""
    if [ -e $INC/request_closest_sql ]; then
      UL1=",ul1=30"
    fi
    $THIS/../trav  -step 1  $INC/search "$QSUB_5$UL1  -N j%n  %Q$THIS/distTree_inc_search_init.sh $INC %f%Q > /dev/null" 
    $THIS/../qstat_wait.sh 2000 1
  fi
  $THIS/distTree_inc_new_log.sh $INC $GRID  
fi


ITER=0
ITER_MAX=`echo $OBJS | awk '{printf "%d", log($1)+3};'`
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
  if [ $GRID_MIN -le 200 ]; then  
    WAIT=7200
  fi
  if [ $GRID == 1 ]; then
    $THIS/../qstat_wait.sh $WAIT 0
  fi
  
  $THIS/distTree_inc_new_log.sh $INC $GRID  

  $THIS/../trav  -step 1  -threads 15  $INC/search "$THIS/distTree_inc_search2bad.sh $INC %f"

  echo "Processing new objects"
  $THIS/distTree_new $QC $INC/  -variance $VARIANCE
done


section "leaf, dissim.add -> tree, dissim"

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
	HYBRID="-hybrid_parent_pairs $INC/hybrid_parent_pairs  -delete_hybrids $INC/hybrid.new  -hybridness_min $HYBRIDNESS_MIN  -dissim_boundary $DISSIM_BOUNDARY  -delete_criterion_outliers $INC/outlier-criterion  -criterion_outlier_num_max 1  -delete_deformation_outliers $INC/outlier-deformation  -deformation_outlier_num_max 1"
fi

DELETE=""
if [ -e $INC/outlier-genogroup ]; then
  wc -l $INC/outlier-genogroup
  DELETE="-delete $INC/outlier-genogroup  -check_delete"
fi

REINSERT=""
POS=$(( ${#VER} - 1 ))
if [ "${VER:$POS}" == 0 ]; then  # PAR
  REINSERT="-reinsert"
  warning "REINSERT\nTree quality temporarily lowers due to missing dissimilarities to close tree neighbors"
fi

# Time: O(n log^4(n)) 
$THIS/makeDistTree $QC  -threads 15  -data $INC/  -variance $VARIANCE \
  $DELETE \
  $REINSERT  -optimize  -skip_len  -subgraph_iter_max 2 \
  -noqual \
  $HYBRID \
  -output_tree $INC/tree.new \
  -dissim_request $INC/dissim_request \
  > $INC/hist/makeDistTree.$VER
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
  $INC/objects_in_tree.sh $INC/outlier-deformation null
  $THIS/../trav $INC/outlier-deformation "$INC/outlier2db.sh %f deformation"  
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


section "Additional requests"
$THIS/distTree_inc_request2dissim.sh $INC $INC/dissim_request $INC/dissim.add-req
if [ -s $INC/dissim.add-req ]; then
  grep -vwi 'nan$' $INC/dissim.add-req | grep -vwi 'inf$' >> $INC/dissim || true
fi
rm $INC/dissim.add-req
rm $INC/dissim_request


$THIS/distTree_inc_tree1_quality.sh $INC


section "QC"
$INC/qc.sh go


echo ""
NEW=`$THIS/distTree_inc_new_list.sh $INC | wc -l`
echo "# New objects: $NEW"
if [ $NEW == 0 -a -z "$REINSERT" ]; then
  touch $INC/finished
else
  rm -f $INC/finished
fi


