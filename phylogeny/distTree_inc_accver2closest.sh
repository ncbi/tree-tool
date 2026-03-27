#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print: <closest accver from #3>\tdistance"
  echo "#1: incremental distance tree directory"
  echo "#2: object accver"
  echo "#3: list of object accvers"
  exit 1
fi
INC=$1
ACCVER=$2
TARGET=$3


#set -x


if grep -wq $ACCVER $TARGET; then
  echo -e "$ACCVER\t0"
  exit
fi


META=$INC/../Metadata.tsv
TREE=$INC/tree.released


# QC
OBJ=$( awk -F "\t" '$2 == "'$ACCVER'"' $META | cut -f 1 )
if [ -z "$OBJ" ]; then
  error "$ACCVER is not found in $META"
fi
$THIS/tree2obj_raw.sh $TREE | grep -wq $OBJ || error "$OBJ is not in $TREE"


TMP=$( mktemp )
#comment $TMP 


echo "#accver" >  $TMP.target
cat $TARGET    >> $TMP.target
$THIS/../tsv/tsv_expand.sh $TMP.target $INC/../Metadata.tsv "" &> /dev/null

tail -n +2 $TMP.target | cut -f 2 | sed 's/^/'$OBJ'\t/1' > $TMP.pair
$THIS/statDistTree $TREE  -dist_request $TMP.pair  -qc  -force  -dist_pairs $TMP.dist  -noprogress > $TMP.statDistTree

L=( $( sort -k3g $TMP.dist | head -1 ) ) || true
CLOSEST=${L[1]}
DIST=${L[2]}
CLOSEST_ACCVER=$( awk -F "\t" '$1 == "'$CLOSEST'"' $META | cut -f 2 )

echo -e "$CLOSEST_ACCVER\t$DIST"


rm $TMP*
