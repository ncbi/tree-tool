#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# -ne 12 ]; then
  echo "Compute dissimilarities"
  echo " #1: file with pairs <object1> <object2>"
  echo " #2: directory with: <object hash>/<object>/<object>.hash-{CDS|PRT}"
  echo " #3: directory with a new object data or ''"
  echo " #4: output file: <object1> <object2> <dissimilarity>"
  echo " #5: hashes intersection_min"
  echo " #6: hashes ratio_min"
  echo " #7: dissim_scale file with dissimilarity thresholds: {CDS PRT univ} | {PRT univ}"
  echo " #8: hmm-univ.stat"
  echo " #9: 1 - BLOSUM62, 0 - PAM30"
  echo "#10: raw power of universal proteins dissimilarity (before averaging)"
  echo "#11: coefficient to multiply the combined dissimilarity by"
  echo "#12: log file (delete on success)"
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
UNIV_POWER=${10}
COEFF=${11}
LOG=${12}


#set -x


CDS=1  # d_CDS is used
DISSIM_SCALE_LINES=`cat $DISSIM_SCALE | grep -v '^ *$' | wc -l`
case $DISSIM_SCALE_LINES in
  2)
    CDS=0
    ;;
  3)
    ;;
  *)
    error "Wrong number of rows in $DISSIM_SCALE"
    ;;
esac

# {CDS,PRT}_RAW_MAX, BARRIER
CDS_RAW_MAX="nan"
BARRIER=1
N=0
if [ $CDS == 1 ]; then
  N=$(( $N + 1 ))
  L=(`head -$N $DISSIM_SCALE | tail -1`)
  CDS_RAW_MAX=${L[1]}
  BARRIER=2
fi
#
N=$(( $N + 1 ))
L=(`head -$N $DISSIM_SCALE | tail -1`)
PRT_RAW_MAX=${L[1]}  # Species barrier
#
if [ $CDS == 1 ]; then
  N=$(( $N + 1 ))
  L=(`head -$N $DISSIM_SCALE | tail -1`)
  #
  N=$(( $N + 1 ))
  L=(`head -$N $DISSIM_SCALE | tail -1`)
fi

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
awk '! ($4 == "nan" || $4 > '$PRT_RAW_MAX')' $TMP.req-dissim-PRT | sed 's/$/\tnan/1' > $TMP.req-dissim-univ_0
awk    '$4 == "nan" || $4 > '$PRT_RAW_MAX    $TMP.req-dissim-PRT > $TMP.req-dissim-PRT_1
req2file $TMP.req-dissim-PRT_1 1 "prot-univ" > $TMP.f1
req2file $TMP.req-dissim-PRT_1 2 "prot-univ" > $TMP.f2
paste $TMP.f1 $TMP.f2 > $TMP.req-univ
$THIS/../../dissim/prot_collection2dissim  -log $LOG  -raw_power $UNIV_POWER  $BLOSUM62  $AVERAGE_MODEL  $TMP.req-univ $TMP.univ
cut -f 3 $TMP.univ > $TMP.dissim-univ_1
paste $TMP.req-dissim-PRT_1 $TMP.dissim-univ_1 > $TMP.req-dissim-univ_1
cat $TMP.req-dissim-univ_0 $TMP.req-dissim-univ_1 > $TMP.req-dissim-univ
  # column 5: d_univ

echo "$TMP.combo ..."
if [ $CDS == 1 ]; then
  mv $TMP.req-dissim-univ $TMP.combo
else
  cut -f 3  --complement $TMP.req-dissim-univ > $TMP.combo
fi

$THIS/../../dissim/combine_dissims $TMP.combo $DISSIM_SCALE  -barrier $BARRIER  -coeff $COEFF  -log $LOG  > $OUT


rm -f $LOG  
rm $TMP*
