#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Dependence of dissimilarity on taxonomy rank"
  echo "#1: dissimilarity file: obj1 obj2 dissim rank"
  exit 1
fi
RPT=$1  


TMP=$( mktemp )
#comment $TMP


function rank2uniKernel ()  
{
  local RANK=$1
  section "rank = $RANK"
  cat $TMP.tab | awk '$2 == '$RANK | cut -f 1 > $TMP
  $THIS/conversion/cols2dm.sh $TMP 6 0 > $TMP.dm
  $THIS/uniKernel $TMP V1 > rank-$RANK.uniKernel
}


# Starts with a BOM
tail -n +3 $RPT | awk '{printf "%f\t%d\n", $3, $4};' > $TMP.tab
cut -f 2 $TMP.tab | sort | uniq -c | awk '$1 >= 1000' | awk '{print $2};' | sort -n > $TMP.rank
#cat $TMP.rank

for RANK in $( cat $TMP.rank )
do
  rank2uniKernel $RANK
done



rm $TMP*
