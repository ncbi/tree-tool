#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print an SQL insert command to populate NCBI.Tax"
  echo "#1: Schema name"
  echo "#2: File with rows: tax_id, parent_tax_id, rank_name, tax_name"
  exit 1
fi
SCHEMA=$1
IN=$2


N=`cat $IN | wc -l`
if [ $N == 0 ]; then
  exit 1
fi
N1=$(( $N - 1 ))

echo "insert into $SCHEMA.Tax (id, parent, rank_name, [name]) values"
head -$N1 $IN | awk -F '\t' '{printf "(%d, %d, ~%s~, ~%s~),\n", $1, $2, $3, $4};' | tr '~' "'"
tail   -1 $IN | awk -F '\t' '{printf "(%d, %d, ~%s~, ~%s~);\n", $1, $2, $3, $4};' | tr '~' "'"
echo "go"


