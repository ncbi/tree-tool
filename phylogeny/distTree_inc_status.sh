#!/bin/csh -f

if ($# != 1) then
  echo "Print criterion statistics for an incremental distance tree"
  echo "#1: distance tree data"
  exit 1
endif


set OBJS = `grep -vc '^ *0x' $1/tree`
echo "# Objects: $OBJS"  

set N = `ls $1/new/ | wc -l`
echo "# To add: $N[1]"

echo ""
grep 'absCriterion =' -n $1/old/makeDistTree.* | sed 's|^'$1'/old/makeDistTree\.||1' | sed 's/:18:/ /1' | grep -v ':19:' | sort -n | head
echo "..."
grep 'absCriterion =' -n $1/old/makeDistTree.* | sed 's|^'$1'/old/makeDistTree\.||1' | sed 's/:18:/ /1' | grep -v ':19:' | sort -n | tail

echo ""
head $1.log
echo "..."
tail $1.log

echo ""
wc -l $1/outlier


