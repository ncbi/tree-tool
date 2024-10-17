#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "#1: incremental distance tree directory"
  echo "#2: new DNA alignment"
  echo "#3: output Newick tree"
  exit 1
fi
INC=$1
ALIGN=$2
OUT=$3


TMP=$( mktemp )
comment $TMP
#set -x


tail -n +2 $ALIGN | cut  -f 1  -d ' ' | sort > $TMP.obj
wc -l $TMP.obj
sort -cu $TMP.obj

$THIS/tree2obj.sh $INC/tree > $TMP.prev.obj
sort -cu $TMP.prev.obj

$THIS/../setMinus $TMP.prev.obj $TMP.obj      > $TMP.del
$THIS/../setMinus $TMP.obj      $TMP.prev.obj > $TMP.new
wc -l $TMP.del
wc -l $TMP.new

if [ -s $TMP.del ]; then
  super_section "Deleting objects"
  $THIS/makeDistTree  -input_tree $INC/tree  -delete $TMP.del  -output_tree $INC/tree.1
  mv $INC/tree.1 $INC/tree
fi


if [ -s $TMP.new ]; then
  super_section "Adding objects"
  QUERY=$TMP.query  
  mkdir $QUERY
  rm -f $INC/query
  ln -s $QUERY $INC/query
  $THIS/../dissim/dna_align_service $ALIGN $QUERY  -threads 20  -log $TMP.query.log  -quiet &
  trap "if [ -e $QUERY ]; then touch $QUERY/.quit; fi" EXIT

  $THIS/distTree_inc_refresh_dissim.sh $INC 0 0

  $THIS/distTree_inc_new_cmd.sh $INC "touch" $TMP.new
  $THIS/distTree_inc.sh $INC 1 '' ''
  rm arc_existence.dm.gz leaf_errors.dm.gz  

  touch $QUERY/.quit
  wait
  rm $QUERY/.quit
  rmdir $QUERY  # QC
  rm $INC/query
  if [ -s $TMP.$QUERY.log ]; then
    error "$TMP.$QUERY.log is not empty"
  fi
fi


$THIS/printDistTree $INC/tree  -order  > $OUT


rm -r $TMP*
