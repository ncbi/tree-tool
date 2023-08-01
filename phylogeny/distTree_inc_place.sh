#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find the place of a sequence in a sequence incremental tree"
  echo "Print: <new sequence name> <tree node> <new leaf arc length> <tree node arc length>"
  echo "#1: incremental tree directory"
  echo "#2: new object file"
  echo "#3: output of statistical report for 100 closest objects"
  exit 1
fi
INC=$1  
QUERY=$2
RES=$3


$INC/qc_object.sh $QUERY


TMP=`mktemp`  
comment $TMP  
#set -x


NAME=`basename $QUERY`
TREE=$INC/tree.released

grep '^ *'$NAME':' $TREE > $TMP.grep || true
if [ -s $TMP.grep ]; then
  rm $TMP*
  error "$NAME is already in $TREE"
fi

section "Finding approximately closest objects"
$INC/object2closest.sh $NAME `dirname $QUERY` | grep -vw "^${NAME}$" | sed 's/$/\t'$NAME'/1' > $TMP.request

section "Fitting to the tree"
cp /dev/null $TMP.dissim
echo "FAIL" > $TMP.leaf
VARIANCE=`cat $INC/variance`
while [ -s $TMP.request ]; do
 #printf "."
  wc -l $TMP.request
  $INC/pairs2dissim.sh $TMP.request $QUERY $TMP.dissim-add $TMP.log  &> /dev/null
  rm $TMP.request
  grep -vwi "nan$" $TMP.dissim-add | grep -vwi "inf$" > $TMP.dissim-add1 || true
  if [ ! -s $TMP.dissim-add1 ]; then
    wc -l $TMP.dissim-add
    warning "Incomparable objects"
    break
  fi
  cat $TMP.dissim-add1 >> $TMP.dissim
  $THIS/distTree_new  $TREE  -variance $VARIANCE  -name $NAME  -dissim $TMP.dissim  -request $TMP.request  -leaf $TMP.leaf  -result $TMP.res
done

section "Placement"
cut -f 1,2,3,4 $TMP.leaf

section "Statistical report"
sort -k3,3g $TMP.res > $TMP.sort
head -100 $TMP.sort > $TMP.res  
mkdir $TMP.report.dir
$THIS/../trav $TMP.res "$INC/pair2report.sh $QUERY 1 %1 > $TMP.report.dir/%n"  -threads 15  
$THIS/../trav $TMP.res "cat $TMP.report.dir/%n" -noprogress > $TMP.report
echo -e "#Object\tObserved dissimilarity\tTree distance\tIdentical proteins\tUniversal proteins compared\tAAI,%" > $RES
paste $TMP.res $TMP.report >> $RES


rm -fr $TMP*

