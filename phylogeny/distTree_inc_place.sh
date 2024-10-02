#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find the place of an object in a released incremental tree"
  echo "Output: #3.placement: <object name> <tree node> <new leaf arc length> <tree node arc length>"
  echo "        #3.neighbors.tsv: statistical report for 100 closest objects"
  echo "#1: incremental tree directory"
  echo "#2: directory with object data"
  echo "#3: output file prefix"
  exit 1
fi
INC=$1  
QUERY=$2
RES=$3


$INC/qc_object.sh $QUERY


NAME=$( basename $QUERY )
TREE=$INC/tree.released


#if false; then 
TMP=$( mktemp )
comment $TMP  
#set -x


grep '^ *'$NAME':' $TREE > $TMP.grep || true
if [ -s $TMP.grep ]; then
  rm $TMP*
  error "$NAME is already in $TREE"
fi

section "Finding approximately closest objects"
$THIS/tree2obj.sh $TREE > $TMP.obj
$INC/object2closest.sh $NAME $QUERY $TMP.obj | grep -vw "^${NAME}$" > $TMP.closest || true
sed 's/$/\t'$NAME'/1' $TMP.closest > $TMP.request

section "Fitting to the tree"
# $TMP.{leaf,closest}
cp /dev/null $TMP.dissim
echo "FAIL" > $TMP.leaf
VARIANCE=$( cat $INC/variance )
while [ -s $TMP.request ]; do
 #printf "."
  wc -l $TMP.request
  $INC/pairs2dissim.sh $TMP.request $QUERY/$NAME $TMP.dissim-add $TMP.log  &> /dev/null
  rm $TMP.request
  grep -vwi "nan$" $TMP.dissim-add | grep -vwi "inf$" > $TMP.dissim-add1 || true
  if [ ! -s $TMP.dissim-add1 ]; then
    wc -l $TMP.dissim-add
    warning "Incomparable objects"
    break
  fi
  cat $TMP.dissim-add1 >> $TMP.dissim
  $THIS/distTree_new  $TREE  -variance $VARIANCE  -name $NAME  -dissim $TMP.dissim  -request $TMP.request  -leaf $TMP.leaf  -closest $TMP.closest  -closest_num 100  -qc
    # PAR
done

section "Placement"
cut -f 1,2,3,4 $TMP.leaf > $RES.placement
cat $RES.placement


section "Statistical report"

echo -e "#Object\tTree distance" >  $TMP.closest.tsv
cat $TMP.closest                 >> $TMP.closest.tsv

# $TMP.neighbors.tsv
cut -f 1 $TMP.closest | sed 's/$/\t'$NAME'/1' > $TMP.request
wc -l $TMP.request
$INC/pairs2dissim.sh $TMP.request $QUERY/$NAME $TMP.dissim-closest $TMP.log  &> /dev/null
echo -e "#Object\tObserved dissimilarity" >  $TMP.neighbors.tsv
sed 's/'$NAME'\t//1' $TMP.dissim-closest  >> $TMP.neighbors.tsv

mkdir $TMP.report.dir
$THIS/../trav -tsv $TMP.neighbors.tsv "$INC/pair2report.sh $QUERY 1 %1 > $TMP.report.dir/%1"  -threads 15  
echo -e "#Object\tIdentical proteins\tUniversal proteins compared\tAAI,%"                     >  $TMP.report.tsv
$THIS/../trav -tsv $TMP.neighbors.tsv "cat $TMP.report.dir/%1 | sed 's/^/%1\t/1'" -noprogress >> $TMP.report.tsv

$THIS/../tsv/tsv_expand.sh $TMP.neighbors.tsv $TMP.closest.tsv -left &> /dev/null
$THIS/../tsv/tsv_expand.sh $TMP.neighbors.tsv $TMP.report.tsv  -left &> /dev/null
if [ -e $INC/../Metadata.tsv ]; then
  cut -f 1,2,3,5,6,9 $INC/../Metadata.tsv > $TMP.meta
  $THIS/../tsv/tsv_expand.sh $TMP.neighbors.tsv $TMP.meta -left &> /dev/null
fi  

$THIS/../tsv/tsv_sort.sh $TMP.neighbors.tsv '-k2,2g -k1,1' > $RES.neighbors.tsv


rm -fr $TMP*

