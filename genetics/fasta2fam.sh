#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Select representatives of a family"
  echo "#1: protein FASTA"
  echo "#2: min. representation share (0..1)"
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


# $TMP.fa
$THIS/fasta2len $F > $TMP.len
cut -f 2 $TMP.len | count > $TMP.count
MEAN=$( grep -w "^mean" $TMP.count | cut -f 2 )
  SD=$( grep -w "^SD"   $TMP.count | cut -f 2 )
# PAR
LEN_MIN=$( echo "$MEAN - 2 * $SD" | bc -l | sed 's/\..*$//1' )
NEG=$( echo "$LEN_MIN < 60" | bc -l )
#
if [ $NEG == 1 ]; then
  echo "MEAN = $MEAN; SD = $SD" 
  > $OUT
else
  $THIS/filterFasta $F -aa  -len_min $LEN_MIN > $TMP.fa

  # $TMP.clust
  # PAR
  $THIS/prot_match.sh $TMP.fa $TMP.fa 0.8 0.8 $CORES > $TMP.match
  $THIS/fa2list.sh $TMP.fa | awk '{OFS="\t"; print $1, $1};' >> $TMP.match
  $THIS/../connectPairs $TMP.match $TMP.clust  -pairs -center

  echo -e "#item\trepr" >  $TMP.clust.tsv
  cat $TMP.clust        >> $TMP.clust.tsv
  $THIS/../tsv/tsv_group $TMP.clust.tsv  -by "repr"  -count "count" > $TMP.clust.grp

  N=$( $THIS/fa2list.sh $TMP.fa | wc -l )
  tail -n +2 $TMP.clust.grp | awk '$2 / '$N' >= '$SHARE | cut -f 1 | sort > $TMP.repr
  sort -cu $TMP.repr

  $THIS/filterFasta $TMP.fa -aa  -target $TMP.repr > $OUT
fi


rm $TMP*  
