#!/bin/csh -f

if ($# != 5) then
  echo "50%-subsampling evaluation of a distance tree"
  echo "#1: Input .dm file without .dm in the parent directory"
  echo "#2: Dissimilarity attribute in #1"
  echo "#3: Random number generator seed, identifier of the item"
  echo "#4: Error log directory"
  echo "#5: Tree building program with parameters #1, #2"
  exit 1
endif


set NAME = `basename $1`
if ($?) exit 1


dm2objs ../$NAME > $NAME-$3.all
if ($?) exit 1

set N = `wc -l $NAME-$3.all`
@ N = $N[1] / 2

setRandOrd $NAME-$3.all  -seed $3 | head -$N > $NAME-$3.subset
set S = $?
if ($S != 141 && $S != 0) then
  echo "subset $S" >> $4/$3
  exit 1
endif
rm $NAME-$3.all

dm2subset ../$NAME $NAME-$3.subset > $NAME-$3.dm
if ($?) exit 1
rm $NAME-$3.subset

$5 $NAME-$3 $2 >> $4/$3
if ($?) then
  echo "Building tree" >> $4/$3
  exit 1
endif
rm $NAME-$3.dm
rm -f $NAME-$3.makeDistTree


rm $4/$3
if ($?) exit 1



