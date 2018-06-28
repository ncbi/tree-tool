#!/bin/csh -f

if ($# != 3) then
  echo "Seed sample of a .dm-file "
  echo "#1: Input .dm file without .dm"
  echo "#2: Sample seed"
  echo "#3: Sample size"
  exit 1
endif

set INPUT     = $1
set BASE_SEED = $2
set BASE_SIZE = $3


set BASE = $INPUT-$BASE_SEED

dm2objs $INPUT | sort > $INPUT.list
if ($?) exit 1

setRandOrd $INPUT.list  -seed $BASE_SEED | head -$BASE_SIZE | sort > $BASE.list
rm $INPUT.list

dm2subset $INPUT $BASE.list > $BASE.dm
if ($?) exit 1
rm $BASE.list


