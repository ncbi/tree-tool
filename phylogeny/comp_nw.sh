#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Comparison of 2 trees with the same leaves by distances"
  echo "#1: tree 1 in Newick"
  echo "#2: tree 2 in Newick"
  echo -e "#3: output file with a line: tree1\ttree2\tleaves\t#dist1>dist2\t#dist2>dist1\t\tmean\tmedian\tmax\tRF_A\tRF_B\tRF"
  echo "#4: output file with pairwise distances between leaves for tree 1and tree 2 (will be gzip'ped) | ''"
  exit 1
fi
T1=$1
T2=$2
OUT=$3
PAIR=$4


$THIS/../check_file.sh $T1 1
$THIS/../check_file.sh $T2 1


TMP=$( mktemp )
comment $TMP


$THIS/newick2tree $T1 > $TMP.tree1
$THIS/newick2tree $T2 > $TMP.tree2

$THIS/tree2obj.sh $TMP.tree1 > $TMP.obj1
$THIS/tree2obj.sh $TMP.tree2 > $TMP.obj2
#wc -l $TMP.obj1
#wc -l $TMP.obj2
diff $TMP.obj1 $TMP.obj2
LEAVES=$( cat $TMP.obj1 | wc -l )


D1=0
D2=0
MEAN=0
MED=0
MAX=0
M1=0
M2=0
diff $T1 $T2 > $TMP.diff || true
if [ -s $TMP.diff ]; then
  super_section "Distance computation"
  $THIS/statDistTree $TMP.tree1  -dist_pairs $TMP.pairs1 1> /dev/null
  $THIS/statDistTree $TMP.tree2  -dist_pairs $TMP.pairs2 1> /dev/null


  super_section "Distance comparison"

  # Object names may include "-"
  sed 's/\t/ - /1' $TMP.pairs1 > $TMP.2col1
  sed 's/\t/ - /1' $TMP.pairs2 > $TMP.2col2
  sort  --parallel 10  -S10G  $TMP.2col1 > $TMP.out1
  sort  --parallel 10  -S10G  $TMP.2col2 > $TMP.out2

  # QC
  cut -f 1 $TMP.out1 > $TMP.label1
  cut -f 1 $TMP.out2 > $TMP.label2
  diff $TMP.label1 $TMP.label2

  if [ $PAIR ]; then
    paste $TMP.out1 $TMP.out2 | cut -f 1,2,4 > $PAIR
    gzip $PAIR &
  fi

  cut -f 2 $TMP.out1 > $TMP.dist1
  cut -f 2 $TMP.out2 > $TMP.dist2
  paste $TMP.dist1 $TMP.dist2 > $TMP.dist

  awk '{print $1 - $2};' $TMP.dist > $TMP.dist_diff_raw
  sort  --parallel 10  -S10G  -g  $TMP.dist_diff_raw > $TMP.dist_diff
  N=$( cat $TMP.dist_diff | wc -l )
  M=$( echo "$N / 2" | bc )
  MED=$( head -$M $TMP.dist_diff | tail -1 )
  cat $TMP.dist_diff | count > $TMP.count
  MEAN=$( grep -w "^mean" $TMP.count | cut -f 2 )
  MAX=$( grep -w "^max" $TMP.count | cut -f 2 )
  D1=$( grep -vx "0" $TMP.dist_diff | grep -c  '^-' || true )
  D2=$( grep -vx "0" $TMP.dist_diff | grep -cv '^-' || true )


  super_section "Robinson-Foulds"
  $THIS/compareTrees $TMP.tree1 $TMP.tree2 > $TMP.comp1
  M1=$( grep -cw '^match-' $TMP.comp1 || true ) 
  $THIS/compareTrees $TMP.tree2 $TMP.tree1 > $TMP.comp2
  M2=$( grep -cw '^match-' $TMP.comp2 || true ) 
  #echo "$M1 $M2"
fi


RF=$(( M1 + M2 ))
echo -e "$T1\t$T2\t$LEAVES\t$D1\t$D2\t$MEAN\t$MED\t$MAX\t$M1\t$M2\t$RF" > $OUT


wait
rm $TMP*
