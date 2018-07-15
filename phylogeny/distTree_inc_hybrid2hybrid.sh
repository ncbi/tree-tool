#!/bin/csh -f

if ($# != 4) then
  echo "Print hybridness information of objects close to a hybrid in the format of makeDistTree if they are hybrid"
  echo "#1: Incremental distance tree directory"
  echo "#2: Hybridness information in the format of makeDistTree"
  echo "#3: Output file"
  echo "#4: log file"
  exit 1
endif


set INC = $1
set HYB = ($2)
set OUT = $3
set LOG = $4

if ($#HYB != 5)  exit 1

set OBJ  = $HYB[1]
set P1   = $HYB[3]
set P2   = $HYB[4]
set D_AB = $HYB[5]


set tmp = `mktemp`
#echo $tmp  


while (1) 
	$INC/db_hybrid2objs.sh $OBJ > $tmp.objs
	if ($? == 0) break
	sleep 30
end
sort $tmp.objs > $tmp.objs1
if ($?) exit 1
mv $tmp.objs1 $tmp.objs
if ($?) exit 1


# Remove close objects which are parents of the hybrid
echo $P1 >  $tmp.parents
echo $P2 >> $tmp.parents
sort $tmp.parents > $tmp.parents1
if ($?) exit 1

setMinus $tmp.objs $tmp.parents1 > $tmp.objs1
if ($?) exit 1
mv $tmp.objs1 $tmp.objs
if ($?) exit 1


set OBJS = `cat $tmp.objs`

if ($#OBJS == 0) goto quit
if ($?) exit 1


trav $tmp.objs "echo %f $P1" >  $tmp.request
if ($?) exit 1
trav $tmp.objs "echo %f $P2" >> $tmp.request
if ($?) exit 1

$INC/request2dissim.sh $tmp.request $tmp.dissim $LOG >& /dev/null
if ($?) exit 1
# Lines of $tmp.dissim may not match those of $tmp.request

# Redundant
cat $tmp.dissim | tr '\t' ' ' > $tmp.dissim1
if ($?) exit 1
mv $tmp.dissim1 $tmp.dissim
if ($?) exit 1

set N = `wc -w $tmp.dissim`
if ($?) exit 1
@ M = $#OBJS * 6
if ($N[1] != $M)  exit 1


set HYBRIDNESS_MIN = `cat $INC/hybridness_min`

while ($#OBJS)
  set OBJ = $OBJS[1]
  grep "^$OBJ " $tmp.dissim >  $tmp.obj_dissim
  grep " $OBJ " $tmp.dissim >> $tmp.obj_dissim  
  set L = `cat $tmp.obj_dissim`
  if ($#L != 6)  exit 1
	set D1 = $L[3]
	set D2 = $L[6]  
	# Cf. distTree_inc_new2hybrid.sh
	set Hybridness = `echo "$D_AB / ($D1 + $D2)" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l`
	set IsHybrid = `echo "$Hybridness >= $HYBRIDNESS_MIN" | bc -l`
	if ($IsHybrid) then
	  echo $OBJS[$i] $Hybridness $P1 $P2 $D_AB | tr ' ' '\t' > $OUT
	  if ($?) exit 1
	endif
	shift OBJS
end



quit:
rm -f $tmp*
rm -f $LOG


