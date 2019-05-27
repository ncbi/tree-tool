#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Compute requested dissimilarities"
  echo "#1: Incremental tree directory"
  echo "#2: Pairs of objects"  
  echo "#3: Output file with dissimilarities"
  exit 1
fi
INC=$1
REQ=$2
OUT=$3


N=`cat $REQ | wc -l`
echo "$N $REQ"
GRID_MIN=`cat $INC/grid_min`
if [ $N -le $GRID_MIN ]; then  # PAR
  $INC/request2dissim.sh $REQ $OUT $OUT.log &> /dev/null
else
	mkdir $INC/dr	
	$THIS/../splitList $REQ $GRID_MIN $INC/dr
	
  while [ 1 == 1 ]; do
    rm -rf $INC/dr.out
		mkdir $INC/dr.out
		
		ls $INC/dr > $INC/dr.list
		while [ -s $INC/dr.list ]; do
  		$THIS/../trav -step 1 $INC/dr.list "$QSUB_5 -N j%f %q$INC/request2dissim.sh $INC/dr/%f $INC/dr.out/%f $INC/dr.out/%f.log%q > /dev/null"		
  		$THIS/../qstat_wait.sh 2000 1	
 			set +o errexit	
		  ls $INC/dr.out/*.log  1> $INC/dr.log  2> /dev/null
	  	set -o errexit
      cat $INC/dr.log | sed 's|^'$INC'/dr.out/||1' | sed 's/\.log$//1' > $INC/dr.list
      rm $INC/dr.log
		done
		rm $INC/dr.list
		
		$THIS/../trav $INC/dr.out "cat %d/%f" > $OUT
		wc -l $OUT
		
		N_new=`cat $OUT | wc -l`
		if [ $N -eq $N_new ]; then 
		  break
		fi
		echo "Redo ..."
	done
	
	rm -r $INC/dr.out
	rm -r $INC/dr
fi	

