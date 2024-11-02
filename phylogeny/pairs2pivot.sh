#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print a pivot table where rows and columns are ordered by distance"
  echo "#1: file with numeric attribute, line format: <row_name> <col_name> <num_attr>"
  echo "#2: output .dm-file"
  exit 1
fi
DATA=$1
OUT=$2


TMP=$( mktemp )
#comment $TMP 


function make_list
{
  local IN=$1
  $THIS/../dm/conversion/obj_attr2dm $IN "Real" -decimals 6 -qc > $TMP.dm
  $THIS/../dm/dm2dist $TMP "dist" -standardize -qc > $TMP.dist.dm
  $THIS/makeDistTree  -data $TMP.dist  -dissim_attr "dist"  -variance "linExp"  -optimize  -subgraph_iter_max 5  -output_tree $TMP.tree >> /dev/stderr
  $THIS/tree2obj_raw.sh $TMP.tree
}


section "Ordering rows" 
make_list $DATA > $TMP.rows

section "Ordering columns" 
awk '{printf "%s\t%s\t%f\n", $2, $1, $3};' $DATA > $TMP.rev
make_list $TMP.rev > $TMP.cols

section "Result"
$THIS/../dm/conversion/obj_attr2dm $DATA "Real" -decimals 6 -qc > $TMP.dm
$THIS/../dm/printDataset $TMP  -objs $TMP.rows  -attrs $TMP.cols > $OUT


rm -f $TMP*


