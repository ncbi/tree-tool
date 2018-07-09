#!/bin/csh -f

if ($# != 5) then
  echo "Print hybrid genome"
  echo "#1: Object id"
  echo "#2: Incremental tree directory"
  echo "#3: min. hybridness (>1)"
  echo "#4: Output file"
  echo "#5: log file"
  exit 1
endif


set tmp = `mktemp`


$2/db2hybrid.sh $1 > $tmp.hyb
if ($?) exit 1

if (-z $tmp.hyb) goto quit


tail -2 $tmp.hyb > $tmp.parents
if ($?) exit 1

set N = `wc -l $tmp.parents`
if ($?) exit 1
if ($N[1] != 2)  exit 1

set P1 = `head -1 $tmp.parents`
set P2 = `tail -1 $tmp.parents`

echo $1 $P1 >  $tmp.request
echo $1 $P2 >> $tmp.request

$2/request2dissim.sh $tmp.request $tmp.dissim $5 >& /dev/null
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
set IsHybrid = `echo "$Hybridness >= $3" | bc -l`
if ($IsHybrid) then
  echo $1 $Hybridness $P1 $P2 $D_AB | tr ' ' '\t' > $4
  if ($N[1] != 2)  exit 1
endif


quit:
rm -f $tmp*
rm -f $5


