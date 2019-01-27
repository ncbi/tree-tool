#!/bin/csh -f

if ($# != 4) then
  echo "Build and evaluate by subsampling a distance tree"
  echo "#1: Input .dm file without .dm"
  echo "#2: Dissimilarity attribute in #1"
  echo "#3: Tree building program with parameters #1, #2"
  echo "#4: #1.tree exists (0/1)"
  exit 1
endif


if (! $4) then
	$3 $1 $2
	if ($?) exit 1
endif


echo ""
echo "Subsampling ..."
# Directory for subsampling .tree-files
mkdir $1.trees
if ($?) exit 1
set DIR = `pwd`
cd $1.trees
if ($?) exit 1

# Error file log
mkdir log
if ($?) exit 1

set replicas = 400  # PAR
set SEED = 0
while ($SEED < $replicas) 
  @ SEED = $SEED + 1
  set OUT = $1-$SEED.tree
  if (! -e $OUT || -z $OUT) then
    cp /dev/null log/$SEED
    if ($?) exit 1
    $QSUB_5 -N j$SEED "subsample_item.sh $1 $2 $SEED log $3" > /dev/null
    if ($?) exit 1
  endif
end
qstat_wait.sh 2000 1
if ($?) exit 1

rmdir log
if ($?) exit 1

cd $DIR


sample_report.sh $1 $replicas none
if ($?) exit 1

#sample_report.sh $1 $replicas directed
#if ($?) exit 1

#sample_report.sh $1 $replicas undirected
#if ($?) exit 1

rm -r $1.trees/  
if ($?) exit 1



