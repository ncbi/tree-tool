#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../bash_common.sh
if [ $# -ne 7 ]; then
  echo "Build a distance tree #2/inc/tree using min_edit distance"
  echo "#1: FASTA"
  echo "#2: initial number of sequences"
  echo "#3: #1 is proteins (0/1)"
  echo "#4: DNA strand is known (0/1)"
  echo "#5: global alignment (0/1)"
  echo "#6: phen/ | ''"
  echo "#7: output directory (to be created)"
  exit 1
fi
FASTA=$( realpath $1 )
INIT_NUM=$2
PROT=$3
KNOWN_STRAND=$4
GLOB=$5
PHEN=$( realpath $6 )
DIR=$7


mkdir $DIR
cd $DIR


TMP=$( mktemp )
comment $TMP


AA=""
DBTYPE="nucl"
if [ $PROT == 1 ]; then
  AA="-aa"
  DBTYPE="prot"
fi

$THIS/../genetics/fa2list.sh $FASTA > $TMP.list
sort -R $TMP.list > $TMP.sort
head -$INIT_NUM $TMP.sort | sort > $TMP.init
$THIS/../genetics/filterFasta $FASTA $AA -whole  -target $TMP.init  -qc > $TMP.fa
$THIS/fasta2tree.sh $TMP.fa $PROT $KNOWN_STRAND $GLOB min_edit init

super_section "inc/"
makeblastdb  -in $FASTA  -dbtype $DBTYPE  -out db  -blastdb_version 4  -logfile $TMP.log
mkdir seq
$THIS/../genetics/splitFasta $FASTA seq $AA -whole -qc
$THIS/distTree_inc_init_stnd.sh inc prot "" "" "" ""
$THIS/../dm/conversion/attr2_2pairs init "dissim" -symmetric -qc | grep -vwi "inf" > inc/dissim
rm init.dm
mv init.tree inc/tree
cp inc/tree inc/hist/tree.1
$THIS/distTree_inc_dissim2indiscern.sh inc inc/dissim
$THIS/../setMinus $TMP.list $TMP.init > $TMP.new
$THIS/distTree_inc_new_cmd.sh inc "touch" $TMP.new
ln -s $PHEN inc/phen
$THIS/tree_quality_phen.sh inc/tree '' inc/phen 0 1 init.qual > inc/hist/tree_quality_phen.1
grep ' !$' inc/hist/tree_quality_phen.1

super_section "Incremental tree building"
mkdir release
$THIS/distTree_inc.sh inc 1 1 "" release
$THIS/printDistTree inc/tree.released -order -qc > release/latest/tree.nw


rm $TMP*
