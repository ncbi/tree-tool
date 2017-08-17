#!/bin/csh -f

if ($# != 1) then
  echo "Bootstrap a feature tree"
  echo "#1: #1.list - list of genomes; #1/ - directory with features of the genomes; #1.tree"
  exit 1
endif



#if (0) then 
# Directory for bootstrap .tree-files
mkdir $1.trees
if ($?) exit 1
cd $1.trees
if ($?) exit 1


# Error file log
mkdir log
if ($?) exit 1
#endif  


set SEED = 0
while ($SEED < 100)  # PAR
  @ SEED = $SEED + 1
  set OUT = $1-$SEED.tree
  if (! -e $OUT || -z $OUT) then
    $QSUB -N j$SEED "featureTree_bootstrap_item.sh $1 $SEED log" > /dev/null
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

ls $1.trees/*.tree | sed 's|^.*/||1' | grep -vw maxParsimony > $1-trees.list
trav $1-trees.list "compareTrees $1.tree $1.trees/%f  -type feature" | grep '^match[+-] ' | sort | uniq -c > $1-time.bootstrap
if ($?) exit 1
rm $1-trees.list
bootstrapReport  $1-time.bootstrap  -replicas 100 | sort -k 4 -n -r
if ($?) exit 1

ls $1.trees/*.tree | sed 's|^.*/||1' | grep -w maxParsimony > $1-trees.list
trav $1-trees.list "compareTrees $1-maxParsimony.tree $1.trees/%f  -type feature" | grep '^match[+-] ' | sort | uniq -c > $1-changes.bootstrap
if ($?) exit 1
rm $1-trees.list
bootstrapReport  $1-changes.bootstrap  -replicas 100 | sort -k 4 -n -r
if ($?) exit 1

