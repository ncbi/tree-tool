#!/bin/csh -f

if ($# != 3) then
  echo "Build and subsampling-bootstrap a distance tree"
  echo "#1: Input .dm file without .dm"
  echo "#2: Distance attribute in #1"
  echo "#3: Tree building program with parameters #1, #2"
  exit 1
endif


$3 $1 $2
if ($?) exit 1


echo ""
echo "Bootstrapping ..."
# Directory for bootstrap .tree-files
mkdir $1.trees
if ($?) exit 1
cd $1.trees
if ($?) exit 1

# Error file log
mkdir log
if ($?) exit 1

while (1)
  set replicas = 400  # PAR
  set SEED = 0
  while ($SEED < $replicas) 
    @ SEED = $SEED + 1
    set OUT = $1-$SEED.tree
    if (! -e $OUT || -z $OUT) then
      cp /dev/null log/$SEED
      if ($?) exit 1
      $QSUB -N j$SEED "bootstrap_item.sh $1 $2 $SEED log $3" > /dev/null
      if ($?) exit 1
    endif
  end
  qstat_wait.sh
  
  rmdir log
  if ($? == 0) break
end


cd ..


bootstrap_report.sh $1 $replicas none
if ($?) exit 1

#bootstrap_report.sh $1 $replicas directed
#if ($?) exit 1

#bootstrap_report.sh $1 $replicas undirected
#if ($?) exit 1

rm -r $1.trees/  
if ($?) exit 1



