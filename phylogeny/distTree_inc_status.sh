#!/bin/csh -f

if ($# != 1) then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  exit 1
endif


set tmp = `mktemp`


echo "Version: `cat $1/version`"
echo "Outliers parameter: `cat $1/strong_outliers`"
echo""

set OBJS = `grep -vc '^ *0x' $1/tree`
echo "# Objects: $OBJS"  

set N = `wc -l $1/leaf`
echo "# Being added: $N[1]"

set N = `ls $1/search/ | wc -l`
echo "# Being searched: $N[1]"

set N = `ls $1/outlier/ | wc -l`
set outliers_percent = `echo "scale=2; $N[1] * 100 / ($N[1] + $OBJS)" | bc -l`
echo "# Outliers: $N[1] ($outliers_percent %)"

echo ""
set N = `ls $1/new/ | wc -l`
echo "# New: $N[1]"

echo ""
wc -l $1/dissim

echo ""
grep '^OUTPUT:' -A 1 -n $1/old/makeDistTree.* | sed 's|^'$1'/old/makeDistTree\.||1' | grep -v ':OUTPUT:' | grep -v '^--$' | sed 's/-[0-9]\+-/ /1' | sort -n -k 1 > $tmp
head -5 $tmp
echo "..."
tail -5 $tmp

echo ""
head -5 $1.log
echo "..."
tail -5 $1.log


rm -f $tmp*
