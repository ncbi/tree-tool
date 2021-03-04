#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print the change of quality of the initial tree"
  echo "#1: Incremental distance tree directory"
  exit 1
fi
INC=$1


if [ ! -e $INC/phen ]; then
  echo "$INC/phen/ must exist"
  exit 1
fi


TMP=`mktemp`
echo $TMP > /dev/stderr

VER=`cat $INC/version`
echo $VER > /dev/stderr
$THIS/tree2obj.sh $INC/hist/tree.1 > $TMP.init
LARGE=0
if [ -e $INC/large ]; then
  LARGE=1
fi
$THIS/tree2obj.sh $INC/tree > $TMP.last
$THIS/../setIntersect.sh $TMP.init $TMP.last 0 > $TMP.list

N=1
while [ $N -lt $VER ]
do
  N=$(( $N + 1 ))
  printf "\r%d" $N > /dev/stderr
  if [ -e $INC/hist/tree.$N.gz ]; then
    gunzip $INC/hist/tree.$N.gz -c > $TMP
  elif [ -e $INC/hist/tree.$N ]; then
    cp $INC/hist/tree.$N $TMP
  else
    continue
  fi
  $THIS/tree_quality_phen.sh $TMP $TMP.list $INC/phen $LARGE 0 "" > $TMP.makeFeatureTree 2> /dev/null
  OBJS=`grep '^# Objects: ' $TMP.makeFeatureTree | sed 's/^#.*: //1'`
  RES=`grep ' !' $TMP.makeFeatureTree | sed 's/^#.*: //1' | sed 's/ .*$//1'`
  echo -e "$N\t$OBJS\t$RES"
done
echo "" > /dev/stderr


rm $TMP*
