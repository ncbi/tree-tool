#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find the place of a sequence in a sequence incremental tree"
  echo "Print: <new sequence name> <tree node> <new leaf arc length> <tree node arc length>"
  echo "#1: incremental tree directory"
  echo "#2: sequence file"
  exit 1
fi
INC=$1  
QUERY=$2


NAME=`basename $QUERY`
ID=`head -1 $QUERY | sed 's/^>//1' | cut -f 1 -d ' '`
if [ "$NAME" != "$ID" ]; then
  error "File name '$NAME' must be the same as FASTA header identifier '$ID'"
fi


TMP=`mktemp`  
#echo $TMP  


$THIS/../genetics/dna_closest.sh $QUERY $INC/seq.fa | grep -vw "$NAME" | sed 's/$/ '$NAME'/1' > $TMP.request

cp /dev/null $TMP.dissim
echo "FAIL" > $TMP.leaf
VARIANCE=`cat $INC/variance`
while [ -s $TMP.request ]; do
  $INC/request2dissim.sh $TMP.request $QUERY $TMP.dissim-add $TMP.log &> /dev/null
  rm $TMP.request
  set +o errexit
  grep -vwi nan $TMP.dissim-add > $TMP.dissim-add1
  set -o errexit
  if [ ! -s $TMP.dissim-add1 ]; then
    break
  fi
  cat $TMP.dissim-add1 >> $TMP.dissim
  $THIS/../phylogeny/distTree_new  $INC/tree.released  -variance $VARIANCE  -name $NAME  -dissim $TMP.dissim  -request $TMP.request  -leaf $TMP.leaf
done

cat $TMP.leaf | cut -f 1,2,3,4


rm -fr $TMP*

