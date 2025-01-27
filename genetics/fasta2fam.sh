#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Select representatives of a family"
  echo "#1: protein FASTA"
  echo "#2: min. representatoin share (0..1)"
  echo "#3: cores"
  echo "#4: output file of representative sequences"
  exit 1
fi
F=$1
SHARE=$2
CORES=$3
OUT=$4


TMP=$( mktemp )
comment $TMP


# PAR
$THIS/prot_match.sh $F $F 0.8 0.8 $CORES > $TMP.match
$THIS/fa2list.sh $F | awk '{OFS="\t"; print $1, $1};' >> $TMP.match
$THIS/../connectPairs $TMP.match $TMP.clust  -pairs -center

echo -e "#item\trepr" >  $TMP.clust.tsv
cat $TMP.clust        >> $TMP.clust.tsv
$THIS/../tsv/tsv_group $TMP.clust.tsv  -by "repr"  -count "count" > $TMP.clust.grp

N=$( $THIS/fa2list.sh $F | wc -l )
tail -n +2 $TMP.clust.grp | awk '$2 / '$N' >= '$SHARE | cut -f 1 | sort > $TMP.repr
sort -cu $TMP.repr
#cut -f 2 $TMP.clust | sort -u > $TMP.repr 

$THIS/filterFasta $F -aa  -target $TMP.repr > $OUT
#$THIS/../phylogeny/fasta2tree.sh $TMP.fa 1 1 1 "min_edit" $TMP.fa
#$THIS/../phylogeny/printDistTree $TMP.fa.tree -order > $TMP.fa.nw


rm $TMP*  
