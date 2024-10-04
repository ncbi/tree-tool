#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Comparison of 2 trees with the same leaves by distances"
  echo "#1: Tree 1 in Newick"
  echo "#2: Tree 2 in Newick"
  exit 1
fi
T1=$1
T2=$2


$THIS/../check_file.sh $T1 1
$THIS/../check_file.sh $T2 1


TMP=$( mktemp )
#comment $TMP


$THIS/newick2tree $T1 > $TMP.tree1
$THIS/newick2tree $T2 > $TMP.tree2

$THIS/tree2obj.sh $TMP.tree1 > $TMP.obj1
$THIS/tree2obj.sh $TMP.tree2 > $TMP.obj2
#wc -l $TMP.obj1
#wc -l $TMP.obj2
diff $TMP.obj1 $TMP.obj2
LEAVES=$( cat $TMP.obj1 | wc -l )

$THIS/statDistTree $TMP.tree1  -dist_pairs $TMP.pairs1 1> /dev/null
$THIS/statDistTree $TMP.tree2  -dist_pairs $TMP.pairs2 1> /dev/null

sed 's/\t/-/1' $TMP.pairs1 | sort > $TMP.out1
sed 's/\t/-/1' $TMP.pairs2 | sort > $TMP.out2

# QC
cut -f 1 $TMP.out1 > $TMP.label1
cut -f 1 $TMP.out2 > $TMP.label2
diff $TMP.label1 $TMP.label2

cut -f 2 $TMP.out1 > $TMP.dist1
cut -f 2 $TMP.out2 > $TMP.dist2
paste $TMP.dist1 $TMP.dist2 > $TMP.dist

awk '{print $1 - $2};' $TMP.dist | sort -g > $TMP.diff
N=$( cat $TMP.diff | wc -l )
M=$( echo "$N / 2" | bc )
MED=$( head -$M $TMP.diff | tail -1 )
cat $TMP.diff | count > $TMP.count
MEAN=$( grep -w "^mean" $TMP.count | cut -f 2 )
MAX=$( grep -w "^max" $TMP.count | cut -f 2 )

$THIS/compareTrees $TMP.tree1 $TMP.tree2 > $TMP.comp1
M1=$( grep -cw '^match-' $TMP.comp1 ) || M1=0
$THIS/compareTrees $TMP.tree2 $TMP.tree1 > $TMP.comp2
M2=$( grep -cw '^match-' $TMP.comp2 ) || M2=0
#echo "$M1 $M2"
RF=$(( M1 + M2 ))


echo -e "$T1\t$T2\t$LEAVES\t$MEAN\t$MED\t$MAX\t$M1\t$M2\t$RF"


rm $TMP*
