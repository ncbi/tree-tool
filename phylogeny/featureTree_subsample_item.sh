#!/bin/csh -f

if ($# != 3) then
  echo "Make a feature tree"
  echo "#1: In the parent directory: #1.list - list of genomes; #1/ - directory with features of the genomes"
  echo "#2: Random number generator seed, identifier of the item"
  echo "#3: Error log directory"
  exit 1
endif



set BASE = `basename $1`
if ($?) exit 1

cp /dev/null $3/$2
if ($?) exit 1


set N = `wc -l ../$BASE.list`
@ N = $N[1] / 2

setRandOrd ../$BASE.list $2 | head -$N > $BASE-$2.list
set S = $?
if ($S != 141 && $S != 0) then
  echo "subset $S" >> $3/$2
  exit 1
endif

featureTree.sh $BASE-$2 ../$BASE > $BASE-$2.out
if ($?) then
  cat $BASE-$2.featureTree >> $3/$2
  exit 1
endif
rm $BASE-$2.out
rm $BASE-$2.dm
rm $BASE-$2.list


rm $3/$2
if ($?) exit 1
