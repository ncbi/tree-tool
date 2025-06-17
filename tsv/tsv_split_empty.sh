#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Split #1 depending on whether the columns of #2 are empty"
  echo "#1: tsv-file"
  echo "#2: comma-separated columns of #1"
  echo "#3: output tsv-file with the rows of #1 where all columns of #2 are non-empty"
  echo "#4: output tsv-file which is a complement to #3"
  exit 1
fi
F=$1
COLS=$2
OUT1=$3
OUT2=$4


#set -x


COND=""
L=( $( echo $COLS | tr ',' ' ' ) )
i=0
while [ $i -lt ${#L[@]} ]; do
  if [ "$COND" ]; then
    COND="$COND && "
  fi
  COND="${COND}\$${L[$i]} != \"\""
  i=$(( i + 1 ))
done

$THIS/tsv_awk.sh $F  '('"$COND"')' > $OUT1
$THIS/tsv_awk.sh $F '!('"$COND"')' | cut -f $COLS --complement > $OUT2


