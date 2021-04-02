#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Add missing requests to closest objects"
  echo "#1: incremental distance tree directory"
  exit 1
fi
INC=$1


TMP=`mktemp`
echo $TMP


mkdir $INC/closest
$THIS/../trav 1000 -zero -start 0 "mkdir $INC/closest/%n" -threads 10

$THIS/tree2obj.sh $INC/tree > $TMP.list
wc -l $TMP.list

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

$THIS/../trav $INC/closest "$THIS/../trav -noprogress %Qgrep -vx %f %d/%f | sed 's/%/\t%f/1'%Q" > $TMP.req


#rm $TMP*
