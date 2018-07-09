#!/bin/csh -f

if ($# != 4) then
  echo "Print hybridness information of a genome in the format of makeDistTree if it is hybrid"
  echo "#1: Incremental distance tree directory"
  echo "#2: Object id"
  echo "#3: Output file"
  echo "#4: log file"
  exit 1
endif


set INC = $1
set OBJ = $2
set OUT = $3
set LOG = $4


set tmp = `mktemp`


set HYBRIDNESS_MIN = `cat $INC/hybridness_min`


while (1) 
	$INC/db2hybrid.sh $OBJ > $tmp.hyb
	if ($? == 0) break
	sleep 30
end

if (-z $tmp.hyb) goto quit


tail -2 $tmp.hyb > $tmp.parents
if ($?) exit 1

set N = `wc -l $tmp.parents`
if ($?) exit 1
if ($N[1] != 2)  exit 1

set P1 = `head -1 $tmp.parents`
set P2 = `tail -1 $tmp.parents`

echo $OBJ $P1 >  $tmp.request
echo $OBJ $P2 >> $tmp.request

$INC/request2dissim.sh $tmp.request $tmp.dissim $LOG >& /dev/null
if ($?) exit 1

set N = `wc -l $tmp.dissim`
if ($?) exit 1
if ($N[1] != 2)  exit 1

set D_AB = `head -1 $tmp.hyb`
set L = `head -1 $tmp.dissim`
set D_AH = $L[3]
set L = `tail -1 $tmp.dissim`
set D_BH = $L[3]

# PAR
set Hybridness = `echo "$D_AB / ($D_AH + $D_BH)" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l`
set IsHybrid = `echo "$Hybridness >= $HYBRIDNESS_MIN" | bc -l`
if ($IsHybrid) then
  echo $OBJ $Hybridness $P1 $P2 $D_AB | tr ' ' '\t' > $OUT
  if ($N[1] != 2)  exit 1
endif


quit:
rm -f $tmp*
rm -f $LOG


