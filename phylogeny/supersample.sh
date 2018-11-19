#!/bin/csh -f

if ($# != 6) then
  echo "Build and evaluate by supersampling a distance tree"
  echo "#1: Input .dm file without .dm"
  echo "#2: Dissimilarity attribute in #1"
  echo "#3: Tree building program with parameters #1, #2"
  echo "#4: Base sample seed"
  echo "#5: Base sample size"
  echo "#6: Super-sample size"
  exit 1
endif

set INPUT      = $1
set ATTR       = $2
set PROG       = $3
set BASE_SEED  = $4
set BASE_SIZE  = $5
set SUPER_SIZE = $6


set BASE = $INPUT-$BASE_SEED

if ($SUPER_SIZE < $BASE_SIZE) then
  echo "Super-sample size < base sample size"
  exit 1
endif

dm2objs $INPUT | sort > $INPUT.list
if ($?) exit 1
set N = `wc -l $INPUT.list`
if ($SUPER_SIZE > $N[1]) then
  echo "Super-sample size > # objects"
  exit 1
endif

setRandOrd $INPUT.list  -seed $BASE_SEED  -sigpipe | head -$BASE_SIZE | sort > $BASE.list

setMinus $INPUT.list $BASE.list > $BASE.rest
if ($?) exit 1
rm $INPUT.list

dm2subset $INPUT $BASE.list > $BASE.dm
if ($?) exit 1

echo "distTriangle ..."
distTriangle $BASE $ATTR
if ($?) exit 1

echo ""
echo ""
echo "$PROG ..."
$PROG $BASE $ATTR
if ($?) exit 1
rm $BASE.dm


echo ""
echo ""
echo "Supersampling ..."
set DIR = `pwd`
# Directory for supersampling .tree-files
mkdir $BASE.trees
if ($?) exit 1
cd $BASE.trees
if ($?) exit 1

# Error file log
mkdir log
if ($?) exit 1

set replicas = 400  # PAR 
set SEED = 0
while ($SEED < $replicas) 
  @ SEED = $SEED + 1
  set OUT = $BASE-$SEED.tree
  if (! -e $OUT || -z $OUT) then
    cp /dev/null log/$SEED
    if ($?) exit 1
    $QSUB_5 -N j$SEED "supersample_item.sh $INPUT $ATTR $BASE_SEED $SEED $SUPER_SIZE log $PROG" > /dev/null
    if ($?) exit 1
  endif
end
qstat_wait.sh 1
if ($?) exit 1

rmdir log
if ($?) exit 1

cd $DIR


rm $BASE.list
rm $BASE.rest


sample_report.sh $BASE $replicas none
if ($?) exit 1

#sample_report.sh $BASE $replicas directed
#if ($?) exit 1

#sample_report.sh $BASE $replicas undirected
#if ($?) exit 1

rm -r $BASE.trees/
if ($?) exit 1



