#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Build a distance tree #2/inc/tree using precomputed pairwise dissimilarities"
  echo "#1: pairwise dissimilarities file where each line has format: object1 object2 dissimilarity. Missing dissimilarities are inf"
  echo "#2: output directory (to be created)"
  echo "Requires: object names should have no dashes"
  exit 1
fi
FULL=$1
DIR=$2


mkdir $DIR
cd $DIR


section "Object list"
awk '{print $1};' $FULL >  tmp
awk '{print $2};' $FULL >> tmp
sort -u tmp > list
rm tmp
if grep '-' list; then
  error "Names should not contain dashes"
fi
wc -l list

section "dissim_full"
awk '{if ($1 < $2)  printf "%s\t%s\t%s\n", $1, $2, $3};' $FULL >  tmp
awk '{if ($1 > $2)  printf "%s\t%s\t%s\n", $2, $1, $3};' $FULL >> tmp
sort -u tmp | sed 's/\t/-/1' > dissim_full
rm tmp

$THIS/distTree_inc_init_stnd.sh inc $THIS/inc/dissim_full "" "" "" ""
touch inc/dissim_full

section "inc/dissim"
sort -R list | head -3000 | sort > list.init || true
$THIS/../list2pairs list.init > dissim_request
rm list.init
inc/pairs2dissim.sh dissim_request "" inc/dissim log
if [ -e log ]; then
  error "File 'log' contains errors"
fi
rm dissim_request
$THIS/distTree_inc_dissim2indiscern.sh inc inc/dissim
$THIS/../dm/conversion/pairs2dm inc/dissim 1 "dissim" 6 -distance > data.dm

section "Initial tree (version: 1)"
VARIANCE=`cat inc/variance`
$THIS/makeDistTree  -threads 5  -data data  -dissim_attr "dissim"  -variance $VARIANCE  -optimize  -subgraph_iter_max 10  -output_tree inc/tree  > inc/hist/makeDistTree-complete.1
rm data.dm
cp inc/tree inc/hist/tree.1

section "inc/new/"
$THIS/tree2obj.sh inc/tree > tree.list
$THIS/../setMinus list tree.list > new.list
rm tree.list
$THIS/distTree_inc_new_cmd.sh inc "touch" new.list
rm new.list

super_section "Incremental tree building"
$THIS/distTree_inc.sh inc 1 ""
