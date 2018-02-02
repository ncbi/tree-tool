#!/bin/csh -f

if ($# != 3) then
  echo "Compute requested dissimilarities"
  echo "#1: incremental tree"
  echo "#2: pairs of objects"  
  echo "#3: output file with dissimilarities"
  exit 1
endif


set N = `wc -l $2`
echo "$N[1] $2"


set M = 2000  # PAR
if ($N[1] <= $M) then  # PAR
  $1/request2dissim.sh $2 $3
  if ($?) exit 1
else
	mkdir $1/dr
	if ($?) exit 1
	
	splitList $2 $M $1/dr
	if ($?) exit 1
	
	mkdir $1/dr.out
	if ($?) exit 1
	
	trav -step 1 $1/dr "$QSUB -N j%f %q$1/request2dissim.sh %d/%f $1/dr.out/%f%q > /dev/null" 
	if ($?) exit 1
	
	qstat_wait.sh
	
	trav -step 100 $1/dr.out "cat %d/%f" > $3
	if ($?) exit 1
	
	set N_new = `wc -l $3`
	if ($N[1] != $N_new[1])  exit 1
	
	rm -r $1/dr.out
	rm -r $1/dr
endif	

