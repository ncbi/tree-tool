#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# -ne 11 ]; then
  echo "Compute dissimilarities"
  echo " #1: file with pairs <object1> <object2>"
  echo " #2: directory with: <object hash>/<object>/<object>.hash-{CDS|PRT}"
  echo " #3: directory with a new object data or ''"
  echo " #4: Output file: <object1> <object2> <dissimilarity>"
  echo " #5: hashes intersection_min"
  echo " #6: hashes ratio_min"
  echo " #7: dissim_scale file with dissimilarity thresholds: [CDS] PRT symbet univ"
  echo " #8: hmm-univ.stat"
  echo " #9: 1 - BLOSUM62, 0 - PAM30"
  echo "#10: power for universal proteins dissimilarity"
  echo "#11: log file (delete on success)"
  exit 1
fi
REQ=$1
GENOME_DIR=$2
NEW_DIR="$3"
OUT=$4
HASH_INTERSECTION_MIN=$5
HASH_RATIO_MIN=$6
DISSIM_SCALE=$7
AVERAGE_MODEL=$8
BLOSUM62_ARG=$9
POWER=${10}
LOG=${11}


CDS=1
DISSIM_SCALE_LINES=`cat $DISSIM_SCALE | grep -v '^ *$' | wc -l`
case $DISSIM_SCALE_LINES in
  0|1|2) 
    error "Too few rows in $DISSIM_SCALE"
    ;;
  3)
    CDS=0
    ;;
  4)
    ;;
  *)
    error "Too many rows in $DISSIM_SCALE"
    ;;
esac

# {CDS,PRT}_DISSIM_MAX, UNIV_DISSIM_AVG
CDS_DISSIM_MAX="nan"
N=0
if [ $CDS == 1 ]; then
  N=$(( $N + 1 ))
  L=(`head -$N $DISSIM_SCALE | tail -1`)
  CDS_DISSIM_MAX=${L[1]}
fi
#
N=$(( $N + 1 ))
L=(`head -$N $DISSIM_SCALE | tail -1`)
PRT_DISSIM_MAX=${L[1]}
#
N=$(( $N + 1 ))
L=(`head -$N $DISSIM_SCALE | tail -1`)
SYMBET_DISSIM_MAX_GLOBAL=`echo "${L[0]} * ${L[1]}" | bc -l`
#
N=$(( $N + 1 ))
L=(`head -$N $DISSIM_SCALE | tail -1`)
UNIV_DISSIM_AVG=`echo "$SYMBET_DISSIM_MAX_GLOBAL / ${L[0]}" | bc -l`
#echo $CDS_DISSIM_MAX $PRT_DISSIM_MAX $SYMBET_DISSIM_MAX_GLOBAL $UNIV_DISSIM_AVG

BLOSUM62=""
if [ $BLOSUM62_ARG == 1 ]; then
  BLOSUM62="-blosum62"
fi


TMP=`mktemp`
#echo $TMP 
#set -x


function req2file
{
  REQ_=$1
  COL=$2  # 1|2
  SUF=$3
  #
  cut -f $COL $REQ_ > $TMP.req2file_req
  $THIS/../../file2hash $TMP.req2file_req -file -append  -log $LOG  -noprogress | awk '{printf "'$GENOME_DIR'/%s/%s/%s.'$SUF'\n", $1, $2, $2};' > $TMP.req2file_col
  if [ $NEW_DIR ]; then
    NAME=`basename $NEW_DIR`
    sed 's|^\(.*/'$NAME'\.'$SUF'\)$|'$NEW_DIR/$NAME'.'$SUF'|1' $TMP.req2file_col
  else
    cat $TMP.req2file_col
  fi
}


awk '{printf "%s\t%s\n", $1, $2};' $REQ > $TMP.req

echo "$TMP.req-PRT ..."
req2file $TMP.req 1 "hash-PRT" > $TMP.f1
req2file $TMP.req 2 "hash-PRT" > $TMP.f2
paste $TMP.f1 $TMP.f2 > $TMP.req-PRT

echo "$TMP.dissim-PRT ..."
$THIS/../../dissim/hash_request2dissim $TMP.req-PRT $TMP.PRT  -intersection_min $HASH_INTERSECTION_MIN  -ratio_min $HASH_RATIO_MIN   -log $LOG 
cut -f 3 $TMP.PRT > $TMP.dissim-PRT

echo "$TMP.dissim-CDS ..."
if [ $CDS == 1 ]; then
  req2file $TMP.req 1 "hash-CDS" > $TMP.f1
  req2file $TMP.req 2 "hash-CDS" > $TMP.f2
  paste $TMP.f1 $TMP.f2 > $TMP.req-CDS
  $THIS/../../dissim/hash_request2dissim $TMP.req-CDS $TMP.CDS  -intersection_min $HASH_INTERSECTION_MIN  -ratio_min $HASH_RATIO_MIN   -log $LOG 
  cut -f 3 $TMP.CDS > $TMP.dissim-CDS
else
  sed 's/^.*$/nan/1' $TMP.req-PRT > $TMP.dissim-CDS
fi

paste $TMP.req $TMP.dissim-CDS $TMP.dissim-PRT > $TMP.req-dissim-PRT
#     1,2      3               4

echo "$TMP.req-dissim-univ ..."
awk '! ($4 == "nan" || $4 > '$PRT_DISSIM_MAX')' $TMP.req-dissim-PRT | sed 's/$/\tnan/1' > $TMP.req-dissim-univ_0
awk    '$4 == "nan" || $4 > '$PRT_DISSIM_MAX    $TMP.req-dissim-PRT > $TMP.req-dissim-PRT_1
req2file $TMP.req-dissim-PRT_1 1 "prot-univ" > $TMP.f1
req2file $TMP.req-dissim-PRT_1 2 "prot-univ" > $TMP.f2
paste $TMP.f1 $TMP.f2 > $TMP.req-univ
$THIS/../../dissim/prot_collection2dissim  -log $LOG  -power $POWER  $BLOSUM62  $AVERAGE_MODEL  $TMP.req-univ $TMP.univ
cut -f 3 $TMP.univ > $TMP.dissim-univ_1
paste $TMP.req-dissim-PRT_1 $TMP.dissim-univ_1 > $TMP.req-dissim-univ_1
cat $TMP.req-dissim-univ_0 $TMP.req-dissim-univ_1 > $TMP.req-dissim-univ
  # column 5: d_univ

# Slowest dissimilarity
echo "$TMP.combo_raw ..."
awk '! (($3 == "nan" || $3 > '$CDS_DISSIM_MAX') && ($5 == "nan" || $5 < '$UNIV_DISSIM_AVG'))' $TMP.req-dissim-univ | sed 's/$/\tnan/1' > $TMP.req-dissim-symbet_0
awk    '($3 == "nan" || $3 > '$CDS_DISSIM_MAX') && ($5 == "nan" || $5 < '$UNIV_DISSIM_AVG')'  $TMP.req-dissim-univ > $TMP.req-dissim-univ_1
cut -f 1,2 $TMP.req-dissim-univ_1 > $TMP.req-symbet_1
wc -l $TMP.req-symbet_1
$THIS/../../trav  -step 1  $TMP.req-symbet_1 "$THIS/symbet.sh $GENOME_DIR %f 2> /dev/null"  -log $LOG > $TMP.dissim-symbet
paste $TMP.req-dissim-univ_1 $TMP.dissim-symbet > $TMP.req-dissim-symbet_1
cat $TMP.req-dissim-symbet_0 $TMP.req-dissim-symbet_1 | awk '{OFS="\t"; print $1,$2,$3,$4,$6,$5};' > $TMP.combo_raw

echo "$TMP.combo ..."
if [ $CDS == 1 ]; then
  mv $TMP.combo_raw $TMP.combo
else
  cut -f 3  --complement $TMP.combo_raw > $TMP.combo
fi

$THIS/../../dissim/combine_dissims $TMP.combo $DISSIM_SCALE  -log $LOG > $OUT


rm -f $LOG  
rm $TMP*
