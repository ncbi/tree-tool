#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Create #1/"
  echo "#1: incremental distance tree directory"
  echo "#2: DNA alignment"
  echo "#3: output tree"
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

sort -R $TMP.obj > $TMP.sort
head -3000 $TMP.sort | sort > $TMP.init
wc -l $TMP.init

$THIS/distTree_inc_init_stnd.sh $INC "dna_align" '' '' '' ''
  # remove hybrids ??


QUERY=$TMP.query  
mkdir $QUERY
ln -s $QUERY $INC/query
$THIS/../dissim/dna_align_service $ALIGN $QUERY  -threads 20  -log $TMP.query.log  -quiet &

mkdir $TMP.work
cd $TMP.work

$THIS/distTree_inc_complete.sh $INC $TMP.init 0
rm data.dm

$THIS/../setMinus $TMP.obj $TMP.init > $TMP.new
wc -l $TMP.new
$THIS/distTree_inc_new_cmd.sh $INC "touch" $TMP.new
$THIS/distTree_inc.sh $INC 1 1 '' ''
rm arc_existence.dm.gz leaf_errors.dm.gz

touch $QUERY/.quit
wait
rm $QUERY/.quit
rmdir $QUERY  # QC
if [ -s $TMP.$QUERY.log ]; then
  error "$TMP.$QUERY.log is not empty"
fi
rm $INC/query


$THIS/printDistTree $INC/tree  -order  > $OUT


rm -r $TMP*
