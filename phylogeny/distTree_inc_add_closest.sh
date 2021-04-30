#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Add missing dissimilarities to closest objects"
  echo "Append: #1/dissim"
  echo "#1: incremental distance tree directory"
  echo "Time: O(n log^5 n)"
  echo "To be run if object2closest.sh changes"
  exit 1
fi
INC=$1


if [ -e $INC/dissim.add-req ]; then
  error "File $INC/dissim.add-req exists"
fi


TMP=`mktemp`
echo $TMP


mkdir $INC/closest
$THIS/../trav 1000 -zero -start 0 "mkdir $INC/closest/%n" -threads 10

$THIS/tree2obj.sh $INC/tree > $TMP.list
wc -l $TMP.list

# Time: O(n log^4 n)
UL1=""
if [ -e $INC/object2closest.sql ]; then
  # PAR
  UL1=",ul1=30"
fi
$THIS/../grid_wait.sh 1
$THIS/../trav  -step 1  $TMP.list "$QSUB_5$UL1  -N j%n  %Q$INC/object2closest.sh %f %q%q > $INC/closest/%h/%f%Q > /dev/null" 
WAIT=2000  # PAR
SEARCH_GRID_MIN=`cat $INC/object2closest.grid`  
if [ $SEARCH_GRID_MIN -le 200 ]; then  
  WAIT=7200
fi
$THIS/../qstat_wait.sh $WAIT 1

$THIS/../trav $INC/closest "$THIS/../trav -noprogress %d/%f %Qgrep -vx %Pf %Pd/%Pf | sed 's/$/\t%Pf/1'%Q" > $TMP.req0
rm -r $INC/closest/ &

# Time: O(n log^5 n)
awk      '$1 < $2'                 $TMP.req0 >  $TMP.req1
awk '{if ($1 > $2) print $2, $1};' $TMP.req0 >> $TMP.req1
sort -u $TMP.req1 > $TMP.req2

# Time: O(n log^3 n)
awk '{if ($1 < $2) print $1, $2};' $INC/dissim >  $TMP.pair
awk '{if ($1 > $2) print $2, $1};' $INC/dissim >> $TMP.pair
sort -u $TMP.pair > $TMP.pair1

$THIS/../setMinus $TMP.req2 $TMP.pair1 > $TMP.req

$THIS/distTree_inc_request2dissim.sh $INC $TMP.req $INC/dissim
wc -l $INC/dissim


rm $TMP*
