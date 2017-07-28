#!/bin/csh -f

if ($# != 3) then
  echo "Make a feature tree"
  echo "#1: In the parent directory: #1.list - list of genomes; #1/ - directory with features of the genomes"
  echo "#2: Random number generator seed, identifier of the item"
  echo "#3: Error log directory"
  exit 1
endif



cp /dev/null $3/$2
if ($?) exit 1


set N = `wc -l ../$1.list`
@ N = $N[1] / 2

setRandOrd ../$1.list $2 | head -$N > $1-$2.list
set S = $?
if ($S != 141 && $S != 0) then
  echo "subset $S" >> $3/$2
  exit 1
endif

featureTree.sh $1-$2 ../$1 > $1-$2.out
if ($?) then
  cat $1-$2.featureTree >> $3/$2
  exit 1
endif
rm $1-$2.out
rm $1-$2.dm
rm $1-$2.list


rm $3/$2
if ($?) exit 1
