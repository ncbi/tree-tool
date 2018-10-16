#!/bin/bash
source bash_common.sh
if [ $# -ne 3 ]; then
  echo "Compute requested dissimilarities"
  echo "#1: Incremental tree"
  echo "#2: Pairs of objects"  
  echo "#3: Output file with dissimilarities"
  exit 1
fi


N=`cat $2 | wc -l`
echo "$N $2"
GRID_MIN=`cat $1/grid_min`
if [ $N -le $GRID_MIN ]; then  # PAR
  $1/request2dissim.sh $2 $3 $3.log &> /dev/null
else
	mkdir $1/dr	
	splitList $2 $GRID_MIN $1/dr
	
  while [ 1 == 1 ]; do
    rm -rf $1/dr.out
		mkdir $1/dr.out
		
		trav -step 1 $1/dr "$QSUB_5 -N j%f %q$1/request2dissim.sh %d/%f $1/dr.out/%f $1/dr.out/%f.log%q > /dev/null"		
		qstat_wait.sh 1
		
		trav $1/dr.out "cat %d/%f" > $3
		wc -l $3
		
		N_new=`cat $3 | wc -l`
		if [ $N -eq $N_new ]; then 
		  break
		fi
		echo "Redo ..."
	done
	
	rm -r $1/dr.out
	rm -r $1/dr
fi	

