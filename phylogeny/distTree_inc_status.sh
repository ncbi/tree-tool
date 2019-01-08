#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  exit 1
fi


TMP=`mktemp`


echo "Version: `cat $1/version`"
echo ""

OBJS=`grep -vc '^ *0x' $1/tree`
echo "# Objects in tree: $OBJS"  

ADDED=`cat $1/leaf | wc -l`
if [ $ADDED -gt 0 ]; then
  echo "# Being added: $ADDED"
fi

SEARCH=`ls $1/search/ | wc -l`
if [ $SEARCH -gt 0 ]; then
  echo "# Being searched: $SEARCH"
fi

N=`ls $1/hybrid/ | wc -l`
PERCENT=`echo "scale=2; $N * 100 / ($N + $OBJS)" | bc -l`
if [ $N -gt 0 ]; then
  echo "# Hybrids: $N ($PERCENT %)"
fi

if [ -e $1/outlier-alien ]; then
	N=`cat $1/outlier-alien | wc -l`
  echo "# Alien-outliers: $N"
fi

if [ -e $1/outlier-dissim ]; then
	N=`cat $1/outlier-dissim | wc -l`
  echo "# Dissimilarity outliers: $N"
fi

if [ -e $1/outlier-genogroup ]; then
	N=`cat $1/outlier-genogroup | wc -l`
  echo "# Genogroup outliers: $N"
  echo "# Objects to be in tree: $(( $OBJS - $N ))"
fi

N=`ls $1/new/ | wc -l`
if [ $N -gt 0 ]; then
  echo "# New: $N"
  echo "# To process: $(( $ADDED + $SEARCH + $N ))"
fi


echo ""
N=`cat $1/dissim | wc -l`
PERCENT=`echo "scale=2; 100 * $N / ($OBJS * ($OBJS - 1))" | bc -l`
echo "# Dissimilarities: $N"
echo "Dissimilarities per object: $PERCENT % of maximum"

echo ""
grep '^OUTPUT:' -A 1 -n $1/hist/makeDistTree*.* | grep -v '\-qc\.' | sed 's|^'$1'/hist/makeDistTree[^.]*\.||1' | grep -v ':OUTPUT:' | grep -v '^--$' | sed 's/-[0-9]\+-/ /1' | sort -n -k 1 > $TMP
tail -5 $TMP

echo ""
tail -5 $1/runlog

echo ""
grep ' V !' $1/hist/makeFeatureTree-tree1.* | sed 's|^'$1'/hist/makeFeatureTree-tree1\.||1' | sed 's/:#/ #/1' | sort -k 1 -n > $TMP
tail -5 $TMP

echo ""
set +o errexit
wc -l inc/hist/outlier-genogroup.* 1> $TMP.out 2> /dev/null
grep -v total $TMP.out |sed 's/^\(.*\)\.\([0-9]\+\)$/\2 \1/1' | sort -n  > $TMP
S=$?
set -o errexit
if [ $S == 0 ]; then
  tail -5 $TMP
fi


rm -f $TMP*



