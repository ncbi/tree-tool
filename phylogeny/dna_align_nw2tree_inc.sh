NOT finished !!!
#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "#1: DNA alignment (absolute pathname)"
  echo "#2: previous Newick tree (absolute pathname)"
  echo "#3: output Newick tree (absolute pathname)"
  exit 1
fi
ALIGN=$1
PREV=$2
OUT=$3


TMP=$( mktemp )
comment $TMP
set -x


mkdir $TMP.work
cd $TMP.work


tail -n +2 $ALIGN | cut  -f 1  -d ' ' | sort > $TMP.obj
wc -l $TMP.obj
sort -cu $TMP.obj

INC=$TMP.inc  
$THIS/distTree_inc_init_stnd.sh $INC "dna_align" '' '' '' ''

$THIS/newick2tree $PREV -qc > $INC/tree
$THIS/tree2obj.sh $INC/tree > $TMP.prev.obj
sort -cu $TMP.prev.obj

$THIS/../setMinus $TMP.prev.obj $TMP.obj      > $TMP.del
$THIS/../setMinus $TMP.obj      $TMP.prev.obj > $TMP.new
wc -l $TMP.del
wc -l $TMP.new

if [ -s $TMP.del ]; then
  $THIS/makeDistTree  -input_tree $INC/tree  -delete $TMP.del  -output_tree $INC/tree.1
  mv $INC/tree.1 $INC/tree
fi


if [ -s $TMP.new ]; then
  QUERY=$TMP.query  
  mkdir $QUERY
  ln -s $QUERY $INC/query
  $THIS/../dissim/dna_align_service $ALIGN $QUERY  -threads 20  -log $TMP.query.log  -quiet &

  $THIS/distTree_inc_refresh_dissim.sh $INC 1 0 0

  $THIS/distTree_inc_new_cmd.sh $INC "touch" $TMP.new
  $THIS/distTree_inc.sh $INC 1 0 '' ''
  rm arc_existence.dm.gz leaf_errors.dm.gz  

  touch $QUERY/.quit
  wait
  rm $QUERY/.quit
  rmdir $QUERY
  if [ -s $TMP.$QUERY.log ]; then
    error "$TMP.$QUERY.log is not empty"
  fi
fi


$THIS/printDistTree $INC/tree  -order  > $OUT


#rm $TMP*
