#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Dependence of #2 on the other nominal columns"
  echo "#1: tsv-file for ANOVA. Column names should be unique w.r.t. space-underscore replacement"
  echo "#2: target column name"
  echo "#3: uniKernel window (0 - default)"
  echo "#4: output directory"
  exit 1
fi
TSV=$1  
TARGET="$2"
WINDOW=$3
OUT=$4


TMP=$( mktemp )
#comment $TMP
#set -x


head -1 $TSV | sed 's/^#//1' | tr '\t' '\n' > $TMP.header
mapfile -t H < $TMP.header
T=$( grep -n -x "$TARGET" $TMP.header | sed 's/:.*$//1' )

I=0
while [ $I -lt ${#H[@]} ]; do
  if [ "${H[$I]}" != "$TARGET" ]; then
    ATTR="${H[$I]}"
    I1=$(( I + 1 ))
    tail -n +2 $TSV | cut -f $I1 | sort -u > $TMP.list
    M=$( < $TMP.list wc -l )
    echo "$ATTR: $M"
    mapfile -t L < $TMP.list
    J=0
    while [ $J -lt ${#L[@]} ]; do
      V="${L[$J]}"
      echo "  $V"
      $THIS/../tsv/tsv_awk.sh $TSV '$'$I1' == "'"$V"'"' | tail -n +2 | cut -f $T > $TMP.data
      $THIS/conversion/cols2dm.sh $TMP.data 0 6 0 > $TMP.dm
      $THIS/uniKernel $TMP "V1" -noprogress  -window $WINDOW  > "$OUT/$ATTR-$V.uniKernel"
      J=$(( J + 1 ))
    done
  fi
  I=$(( I + 1 ))
done


rm $TMP*
