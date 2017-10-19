#!/bin/csh -f

if ($# != 1) then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  exit 1
endif


echo "Version: `cat $1/version`"
echo""

set OBJS = `grep -vc '^ *0x' $1/tree`
echo "# Objects: $OBJS"  

set N = `ls $1/new/ | wc -l`
echo "# To add: $N[1]"

set N = `ls $1/search/ | wc -l`
echo "# Being searched: $N[1]"

set N = `wc -l $1/leaf`
echo "# Being added: $N[1]"

wc -l $1/outlier

echo ""
wc -l $1/dissim

echo ""
grep 'absCriterion =' -n $1/old/makeDistTree.* | sed 's|^'$1'/old/makeDistTree\.||1' | sed 's/:18:/ /1' | grep -v ':19:' | sort -n | head
echo "..."
grep 'absCriterion =' -n $1/old/makeDistTree.* | sed 's|^'$1'/old/makeDistTree\.||1' | sed 's/:18:/ /1' | grep -v ':19:' | sort -n | tail

echo ""
head $1.log
echo "..."
tail $1.log



