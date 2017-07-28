#!/bin/csh -f

if ($# != 2) then
  echo "Build and bootstrap a distance tree"
  echo "#1: Input .dm file without .dm; #1.tree must exist"
  echo "#2: Distance attribute in #1"
  exit 1
endif


distTree.sh $1 $2
if ($?) exit 1


echo ""
echo "Bootstrapping..."
#if (0) then  # ??
# Directory for bootstrap .tree-files
mkdir $1.trees
if ($?) exit 1
cd $1.trees
if ($?) exit 1


# Error file log
mkdir log
if ($?) exit 1
#endif  # ??


set SEED = 0
while ($SEED < 100)  # PAR
  @ SEED = $SEED + 1
  set OUT = $1-$SEED.tree
  if (! -e $OUT || -z $OUT) then
    $QSUB -N j$SEED "distTree_bootstrap_item.sh $1 $2 $SEED log" > /dev/null
  endif
end
while (1)
  sleep 15  # PAR
  set Q = `qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | head -1 | wc -l`
  if ($Q[1] == 0)  break
end

rmdir log
if ($?) exit 1


cd ..

bootstrap.sh $1
if ($?) exit 1


