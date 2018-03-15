#!/bin/csh -f

if ($# != 7) then
  echo "#1: Input .dm file without .dm in the parent directory, #1-#4.list, #1-#4.rest"
  echo "#2: Dissimilarity attribute in #1"
  echo "#3: Base sample seed"
  echo "#4: Sample seed"
  echo "#5: Sample size"
  echo "#6: Error log directory"
  echo "#7: Tree building program with parameters #1, #2"
  exit 1
endif

set INPUT     = `basename $1`
set ATTR      = $2
set BASE_SEED = $3
set SEED      = $4
set SIZE      = $5
set LOG       = $6
set PROG      = $7


set BASE = $INPUT-$BASE_SEED

set BASE_SIZE = `wc -l ../$BASE.list`
if ($SIZE < $BASE_SIZE[1]) then
  echo "Sample size < base sample size"
  exit 1
endif
@ N = $SIZE - $BASE_SIZE[1]

cp ../$BASE.list $BASE-$SEED.list
if ($?) exit 1
setRandOrd ../$BASE.rest $SEED | head -$N >> $BASE-$SEED.list
set S = $?
if ($S != 141 && $S != 0) then
  echo "subset $S" >> $LOG/$SEED
  exit 1
endif

dm2subset ../$INPUT $BASE-$SEED.list > $BASE-$SEED.dm
if ($?) exit 1
rm $BASE-$SEED.list

$PROG $BASE-$SEED $ATTR >> $LOG/$SEED
if ($?) then
  echo "Building tree" >> $LOG/$SEED
  exit 1
endif
rm $BASE-$SEED.dm
rm -f $BASE-$SEED.makeDistTree


rm $LOG/$SEED
if ($?) exit 1



