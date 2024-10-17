#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  echo "#2: number of latest versions to report"
  exit 1
fi
INC=$1
HIST=$2


echo "Version: $( cat $INC/version )"

set +o errexit
OBJS=$( grep -vc '^ *0x' $INC/tree )
set -o errexit
section "# Objects in tree: $OBJS"  

if [ $OBJS == 0 ]; then
  exit
fi


TMP=$( mktemp )


ADDED=$( cat $INC/leaf | wc -l )
if [ $ADDED -gt 0 ]; then
  echo "# Being added: $ADDED"
fi

SEARCH=$( ls $INC/search/ | wc -l )
if [ $SEARCH -gt 0 ]; then
  echo "# Being searched: $SEARCH"
fi


N=$( $THIS/distTree_inc_new_list.sh $INC | wc -l )
if [ $N -gt 0 ]; then
  section "# New: $N"
  echo "# To process: $(( $ADDED + $SEARCH + $N ))"
fi

if [ -e $INC/outlier-genogroup ]; then
	N=$( cat $INC/outlier-genogroup | wc -l )
  section "# Genogroup outliers: $N"
  echo "# Objects to be in tree: $(( $OBJS - $N ))"
fi


echo ""
N=( $( wc -l $INC/dissim ) )
PERCENT=$( echo "scale=2; 100 * ${N[0]} / ($OBJS * ($OBJS - 1) / 2)" | bc -l )  # May be > 100% due to outliers which have dissimilarities
echo "# Dissimilarities: ${N[0]}"
echo "Dissimilarities per object: $PERCENT % of maximum"

echo ""
grep '^OUTPUT:' -A 1 -n $INC/hist/makeDistTree*.* | grep -v '\-qc\.' | sed 's|^'$INC'/hist/makeDistTree[^.]*\.||1' | grep -v ':OUTPUT:' | grep -v '^--$' | sed 's/-[0-9]\+-/ /1' | sort -n -k 1 > $TMP
tail -$HIST $TMP

echo ""
tail -$HIST $INC/runlog

set +o errexit
grep ' V !$' $INC/hist/makeFeatureTree-tree1.* 1> $TMP.grep 2> /dev/null
set -o errexit
if [ -s $TMP.grep ]; then
  sed 's|^'$INC'/hist/makeFeatureTree-tree1\.||1' $TMP.grep | sed 's/:#/ #/1' | sort -k 1 -n > $TMP 
  echo ""
  tail -$HIST $TMP
fi


function report_count
{
  local NAME="$1"
  set +o errexit
  wc -l $INC/hist/$NAME.* 1> $TMP.out 2> /dev/null
  grep -v "total" $TMP.out | sed 's/^\(.*\)\.\([0-9]\+\)$/\2 \1/1' | sort -n  > $TMP
  S=$?
  set -o errexit
  if [ $S == 0 ]; then
    echo ""
    tail -$HIST $TMP
  fi
}


#report_count "hybrid"  
report_count "hybrid-indiscern"
#report_count "unhybrid"
report_count "outlier-genogroup"


rm -f $TMP*



