#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  echo "#2: number of latest versions to report"
  exit 1
fi
INC=$1
HIST=$2


echo "Version: `cat $INC/version`"
echo ""

set +o errexit
OBJS=`grep -vc '^ *0x' $INC/tree`
set -o errexit
echo "# Objects in tree: $OBJS"  

if [ $OBJS == 0 ]; then
  exit
fi


TMP=`mktemp`


ADDED=`cat $INC/leaf | wc -l`
if [ $ADDED -gt 0 ]; then
  echo "# Being added: $ADDED"
fi

SEARCH=`ls $INC/search/ | wc -l`
if [ $SEARCH -gt 0 ]; then
  echo "# Being searched: $SEARCH"
fi


N=`ls $INC/new/ | wc -l`
if [ $N -gt 0 ]; then
  echo ""
  echo "# New: $N"
  echo "# To process: $(( $ADDED + $SEARCH + $N ))"
fi

if [ -e $INC/outlier-genogroup ]; then
	N=`cat $INC/outlier-genogroup | wc -l`
  echo ""
  echo "# Genogroup outliers: $N"
  echo "# Objects to be in tree: $(( $OBJS - $N ))"
fi


echo ""
N=`cat $INC/dissim | wc -l`
PERCENT=`echo "scale=2; 100 * $N / ($OBJS * ($OBJS - 1) / 2)" | bc -l`  # May be > 100% due to outliers which have dissimilarities
echo "# Dissimilarities: $N"
echo "Dissimilarities per object: $PERCENT % of maximum"

echo ""
grep '^OUTPUT:' -A 1 -n $INC/hist/makeDistTree*.* | grep -v '\-qc\.' | sed 's|^'$INC'/hist/makeDistTree[^.]*\.||1' | grep -v ':OUTPUT:' | grep -v '^--$' | sed 's/-[0-9]\+-/ /1' | sort -n -k 1 > $TMP
tail -$HIST $TMP

echo ""
tail -$HIST $INC/runlog

set +o errexit
grep ' V !' $INC/hist/makeFeatureTree-tree1.* 1> $TMP.grep 2> /dev/null
set -o errexit
if [ -s $TMP.grep ]; then
  sed 's|^'$INC'/hist/makeFeatureTree-tree1\.||1' $TMP.grep | sed 's/:#/ #/1' | sort -k 1 -n > $TMP 
  echo ""
  tail -$HIST $TMP
fi

set +o errexit
wc -l $INC/hist/hybrid.* 1> $TMP.out 2> /dev/null
grep -v total $TMP.out | sed 's/^\(.*\)\.\([0-9]\+\)$/\2 \1/1' | sort -n  > $TMP
S=$?
set -o errexit
if [ $S == 0 ]; then
  echo ""
  tail -$HIST $TMP
fi

if [ 0 == 1 ]; then
  set +o errexit
  wc -l $INC/hist/unhybrid.* 1> $TMP.out 2> /dev/null
  grep -v total $TMP.out | sed 's/^\(.*\)\.\([0-9]\+\)$/\2 \1/1' | sort -n  > $TMP
  S=$?
  set -o errexit
  if [ $S == 0 ]; then
    echo ""
    tail -$HIST $TMP
  fi
fi

set +o errexit
wc -l $INC/hist/outlier-genogroup.* 1> $TMP.out 2> /dev/null
grep -v total $TMP.out | sed 's/^\(.*\)\.\([0-9]\+\)$/\2 \1/1' | sort -n  > $TMP
S=$?
set -o errexit
if [ $S == 0 ]; then
  echo ""
  tail -$HIST $TMP
fi


rm -f $TMP*



