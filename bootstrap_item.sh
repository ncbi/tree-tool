#!/bin/csh -f

if ($# != 5) then
  echo "50%-Bootstrap a distance tree"
  echo "#1: Input .dm file without .dm in the parent directory"
  echo "#2: Distance attribute in #1"
  echo "#3: Random number generator seed, identifier of the item"
  echo "#4: Error log directory"
  echo "#5: Tree building program with parameters #1, #2"
  exit 1
endif



cp /dev/null $4/$3
if ($?) exit 1


dm2objs ../$1 > $1-$3.all
if ($?) exit 1

set N = `wc -l $1-$3.all`
@ N = $N[1] / 2

setRandOrd $1-$3.all $3 | head -$N > $1-$3.subset
set S = $?
if ($S != 141 && $S != 0) then
  echo "subset $S" >> $4/$3
  exit 1
endif
rm $1-$3.all

dm2subset ../$1 $1-$3.subset > $1-$3.dm
if ($?) exit 1
rm $1-$3.subset

$5 $1-$3 $2 >> $4/$3
if ($?) then
  echo "Building tree" >> $4/$3
  exit 1
endif
rm $1-$3.dm
rm -f $1-$3.makeDistTree


rm $4/$3
if ($?) exit 1



