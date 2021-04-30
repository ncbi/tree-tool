#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find the place of a sequence in a sequence incremental tree"
  echo "Print: <new sequence name> <tree node> <new leaf arc length> <tree node arc length>"
  echo "#1: incremental tree directory"
  echo "#2: new object file or directory"
  exit 1
fi
INC=$1  
QUERY=$2


$INC/qc_object.sh $QUERY


TMP=`mktemp`  
#echo $TMP  
#set -x

NAME=`basename $QUERY`

grep '^ *'$NAME':' $INC/tree > $TMP.grep || true
if [ -s $TMP.grep ]; then
  rm $TMP*
  error "$NAME is already in $INC/tree"
fi

section "Finding closest objects"
$INC/object2closest.sh $NAME `dirname $QUERY` | grep -vw "^${NAME}$" | sed 's/$/\t'$NAME'/1' > $TMP.request

section "Fitting to the tree"
cp /dev/null $TMP.dissim
echo "FAIL" > $TMP.leaf
VARIANCE=`cat $INC/variance`
while [ -s $TMP.request ]; do
  printf "."
  $INC/pairs2dissim.sh $TMP.request $QUERY $TMP.dissim-add $TMP.log  &> /dev/null
  rm $TMP.request
  grep -vwi "nan$" $TMP.dissim-add | grep -vwi "inf$" > $TMP.dissim-add1 || true
  if [ ! -s $TMP.dissim-add1 ]; then
    wc -l $TMP.dissim-add
    warning "Incomparable objects"
    break
  fi
  cat $TMP.dissim-add1 >> $TMP.dissim
  $THIS/distTree_new  $INC/tree  -variance $VARIANCE  -name $NAME  -dissim $TMP.dissim  -request $TMP.request  -leaf $TMP.leaf
done
echo ""
echo ""

cat $TMP.leaf | cut -f 1,2,3,4


rm -fr $TMP*

