#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 12 ]; then
  echo "Compute dissimilarities"
  echo " #1: file with pairs <Object1> <Object2>"
  echo " #2: directory with: <Object hash>/<Object>/<Object>.hash-{CDS|PRT}"
  echo " #3: directory with a new object data or ''"
  echo " #4: Output file"
  echo " #5: CDS is used (0/1)"
  echo " #6: hashes intersection_min"
  echo " #7: hashes ratio_min"
  echo " #8: dissim_scale"
  echo " #9: hmm-univ.stat"
  echo "#10: 1 - BLOSUM62, 0 - PAM30"
  echo "#11: power for universal proteins dissimilarity"
  echo "#12: log file (delete on success)"
  exit 1
fi
IN=$1
GENOME_DIR=$2
NEW_DIR="$3"
OUT=$4
CDS=$5
HASH_INTERSECTION_MIN=$6
HASH_RATIO_MIN=$7
DISSIM_SCALE=$8
AVERAGE_MODEL=$9
BLOSUM62_ARG=${10}
POWER=${11}
LOG=${12}


DISSIM_SCALE_LINES=2
if [ $CDS == 1 ]; then
  DISSIM_SCALE_LINES=3
fi
N=`cat $DISSIM_SCALE | wc -l`
if [ $N -ne $DISSIM_SCALE_LINES ]; then
  echo "$DISSIM_SCALE must have $DISSIM_SCALE_LINES lines" >> $LOG
  exit 1
fi

PRT_LINE=1
if [ $CDS == 1 ]; then
  PRT_LINE=2
fi
L=(`head -$PRT_LINE $DISSIM_SCALE | tail -1`)
DISSIM_MAX=${L[2]}

BLOSUM62=""
if [ $BLOSUM62_ARG == 1 ]; then
  BLOSUM62="-blosum62"
fi


TMP=`mktemp`
#echo $TMP 
#set -x


# $TMP.req
cat $IN | awk '{print $1};' > $TMP.1
cat $IN | awk '{print $2};' > $TMP.2
$THIS/../trav $TMP.1 "echo $GENOME_DIR/%h/%f/%f.hash-PRT" > $TMP.f1
$THIS/../trav $TMP.2 "echo $GENOME_DIR/%h/%f/%f.hash-PRT" > $TMP.f2
paste $TMP.f1 $TMP.f2 | tr '\t' ' ' > $TMP.req
if [ $NEW_DIR ]; then
  NAME=`basename $NEW_DIR`
  sed 's|^\(.*/'$NAME'\.hash-PRT\) |'$NEW_DIR/$NAME'.hash-PRT |1' $TMP.req | sed 's| \(.*/'$NAME'\.hash-PRT\)$| '$NEW_DIR/$NAME'.hash-PRT|1' > $TMP.req-mod
  mv $TMP.req-mod $TMP.req
fi

# $TMP.prt{0|1}
$THIS/hash_request2dissim $TMP.req $TMP.prt  -intersection_min $HASH_INTERSECTION_MIN  -ratio_min $HASH_RATIO_MIN   -log $LOG 
set +o errexit
grep -w  "nan" $TMP.prt                           >  $TMP.prt1
grep -wv "nan" $TMP.prt | awk '$3 >  '$DISSIM_MAX >> $TMP.prt1
grep -wv "nan" $TMP.prt | awk '$3 <= '$DISSIM_MAX >  $TMP.prt0
set -o errexit

# $TMP.prt0 -> $TMP.{cds,univ}0
cat $TMP.prt0 | tr '\t' ' ' | sed 's/ [^ ]\+$/ /1' | sed 's/\.hash-PRT /.hash-CDS /g' > $TMP.req0
if [ $CDS == 1 ]; then
  $THIS/hash_request2dissim $TMP.req0 $TMP.cds0  -intersection_min $HASH_INTERSECTION_MIN  -ratio_min $HASH_RATIO_MIN   -log $LOG 
fi
sed 's/$/ nan/1' $TMP.req0 > $TMP.univ0

# $TMP.prt1 -> $TMP.{cds,univ}1
cat $TMP.prt1 | tr '\t' ' ' | sed 's/ [^ ]\+$/ /1' | sed 's/\.hash-PRT /.prot-univ /g' > $TMP.req1
if [ $CDS == 1 ]; then
  sed 's/$/ nan/1' $TMP.req1 > $TMP.cds1
fi
$THIS/prots_pair2dissim  -log $LOG  -power $POWER  $BLOSUM62  $AVERAGE_MODEL $TMP.req1 $TMP.univ1

# $TMP.combo-raw
if [ $CDS == 1 ]; then
  cat $TMP.cds0 $TMP.cds1 $TMP.prt0 $TMP.prt1 $TMP.univ0 $TMP.univ1 > $TMP.combo-raw
else
  cat                     $TMP.prt0 $TMP.prt1 $TMP.univ0 $TMP.univ1 > $TMP.combo-raw
fi

# $TMP.combo
cat $TMP.combo-raw | tr '\t' ' ' | sed 's|'$GENOME_DIR'/[^/]\+/[^/]\+/||g' | sed 's/\.hash-... / /g' | sed 's/\.prot-univ / /g' > $TMP.combo
if [ $NEW_DIR ]; then
  sed 's|^'$NEW_DIR'/||1' $TMP.combo | sed 's| '$NEW_DIR'/| |1' > $TMP.combo-mod
  mv $TMP.combo-mod $TMP.combo
fi

$THIS/combine_dissims $TMP.combo $DISSIM_SCALE  -log $LOG > $OUT


rm -f $LOG  
rm $TMP*
