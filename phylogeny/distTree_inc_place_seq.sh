#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find the place of a sequence in a sequecne incremewntal tree"
  echo "#1: sequence file"
  echo "#2: incremental tree directory"
  exit 1
fi
QUERY=$1
INC=$2  


TMP=`mktemp`  
#echo $TMP  


NAME=`basename $QUERY`
$THIS/../genetics/dna_closest.sh $QUERY $INC/seq.fa | sed 's/$/ '$NAME'/1' > $TMP.request

cp /dev/null $TMP.dissim
echo "FAIL" > $TMP.leaf
VARIANCE=`cat $INC/variance`
while [ -s $TMP.request ]; do
  $INC/request2dissim.sh $TMP.request $QUERY $TMP.dissim-add $TMP.log &> /dev/null
  rm $TMP.request
  cat $TMP.dissim-add >> $TMP.dissim
  $THIS/../phylogeny/distTree_new  $INC/  -variance $VARIANCE  -name $NAME  -dissim $TMP.dissim  -request $TMP.request  -leaf $TMP.leaf
done

L=(`cat $TMP.leaf`)
echo "${L[0]} has arc of length ${L[2]} joining above ${L[1]} by ${L[3]}"


rm -fr $TMP*

