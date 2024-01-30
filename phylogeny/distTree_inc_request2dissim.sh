#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
source $THIS/../qsub_env.sh
if [ $# -ne 3 ]; then
  echo "Compute requested dissimilarities"
  echo "#1: incremental tree directory"
  echo "#2: pairs of objects"  
  echo "#3: output file with dissimilarities"
  exit 1
fi
INC=$1
REQ=$2
OUT=$3


N=`cat $REQ | wc -l`
echo "$N $REQ"


# $INC/dr.res
GRID_MIN=`cat $INC/pairs2dissim.grid`
if [ -e $INC/dissim_full -o $N -le $GRID_MIN ]; then 
  $INC/pairs2dissim.sh $REQ "" $INC/dr.res $OUT.log > /dev/null
  if [ -e $OUT.log ]; then
    if [ -s $OUT.log ]; then
      head $OUT.log
      exit 1
    fi
    rm $OUT.log
  fi
else
	mkdir $INC/dr	
	TMP=`mktemp`
	sort -R $REQ > $TMP
	$THIS/../splitList $TMP $GRID_MIN $INC/dr
	rm $TMP

  # $INC/dr.{out/,res}
  while true; do
    rm -rf $INC/dr.out
		mkdir $INC/dr.out
		
		ls $INC/dr > $INC/dr.list
		while [ -s $INC/dr.list ]; do
		  wc -l $INC/dr.list
		  mkdir $INC/dr.log
		  $THIS/../grid_wait.sh 1
  		$THIS/../trav  -step 1  $INC/dr.list "$QSUB_5 -N j%f %q$INC/pairs2dissim.sh $INC/dr/%f %Q%Q $INC/dr.out/%f $INC/dr.log/%f%q > /dev/null"		
  		$THIS/../qstat_wait.sh 7200 1	 
		  ls $INC/dr.log > $INC/dr.bad
		  ls $INC/dr.out > $INC/dr.good
	  	$THIS/../setMinus $INC/dr.list $INC/dr.good >> $INC/dr.bad
	  	rm $INC/dr.good
	  	sort -u $INC/dr.bad > $INC/dr.list
	  	rm $INC/dr.bad
	  	rm -r $INC/dr.log
		done
		rm $INC/dr.list
		
		$THIS/../trav $INC/dr.out "cat %d/%f" > $INC/dr.res

		N_new=`cat $INC/dr.res | wc -l`
		echo "$N_new $INC/dr.res"
		if [ $N -eq $N_new ]; then 
		  break
		fi
		warning "Redo"
	done
	
	rm -r $INC/dr.out &
	rm -r $INC/dr &
fi	


grep -vwi 'nan$' $INC/dr.res | grep -vwi 'inf$' > $OUT || true
rm $INC/dr.res
wc -l $OUT


wait
