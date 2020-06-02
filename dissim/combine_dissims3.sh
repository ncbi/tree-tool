#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 10 ]; then
  echo "Compute dissimilarities"
  echo "#1: file with pairs <Object1> <Object2>"
  echo "#2: Output file"
  echo "#3: directory with: <Object hash>/<Object>/<Object>.hash-{CDS|PRT}"
  echo "#4: hashes intersection_min"
  echo "#5: hashes ratio_min"
  echo "#6: dissim_scale"
  echo "#7: hmm-univ.stat"
  echo "#8: 1 - BLOSUM62, 0 - PAM30"
  echo "#9: power for universal proteins dissimilarity"
  echo "#10: log file (delete on success)"
  exit 1
fi
IN=$1
OUT=$2
GENOME_DIR=$3
HASH_INTERSECTION_MIN=$4
HASH_RATIO_MIN=$5
DISSIM_SCALE=$6
AVERAGE_MODEL=$7
BLOSUM62_ARG=$8
POWER=$9
LOG=${10}


N=`cat $DISSIM_SCALE | wc -l`
if [ $N -ne 3 ]; then
  echo "$DISSIM_SCALE must have 3 lines" >> $LOG
  exit 1
fi

L=(`head -2 $DISSIM_SCALE | tail -1`)
DISSIM_MAX=${L[2]}

BLOSUM62=""
if [ $BLOSUM62_ARG == 1 ]; then
  BLOSUM62="-blosum62"
fi


TMP=`mktemp`
#echo $TMP 


cat $IN | awk '{print $1};' > $TMP.1
cat $IN | awk '{print $2};' > $TMP.2
$THIS/../trav $TMP.1 "echo %n $GENOME_DIR/%h/%f/%f.hash-PRT" > $TMP.f1
$THIS/../trav $TMP.2 "echo %n $GENOME_DIR/%h/%f/%f.hash-PRT" > $TMP.f2
join  -1 1  -2 1  $TMP.f1 $TMP.f2 | cut -d ' ' -f 2,3 > $TMP.req

# $TMP.prt{0|1}
$THIS/hash_request2dissim $TMP.req $TMP.prt  -intersection_min $HASH_INTERSECTION_MIN  -ratio_min $HASH_RATIO_MIN   -log $LOG 
set +o errexit
grep -w  nan $TMP.prt                           >  $TMP.prt1
grep -wv nan $TMP.prt | awk '$3 >  '$DISSIM_MAX >> $TMP.prt1
grep -wv nan $TMP.prt | awk '$3 <= '$DISSIM_MAX >  $TMP.prt0
set -o errexit

# $TMP.prt0 -> $TMP.{cds,univ}0
cat $TMP.prt0 | tr '\t' ' ' | sed 's/ [^ ]\+$/ /1' | sed 's/\.hash-PRT /.hash-CDS /g' > $TMP.req0
$THIS/hash_request2dissim $TMP.req0 $TMP.cds0  -intersection_min $HASH_INTERSECTION_MIN  -ratio_min $HASH_RATIO_MIN   -log $LOG 
sed 's/$/ nan/1' $TMP.req0 > $TMP.univ0

# $TMP.prt1 -> $TMP.{cds,univ}1
cat $TMP.prt1 | tr '\t' ' ' | sed 's/ [^ ]\+$/ /1' | sed 's/\.hash-PRT /.prot-univ /g' > $TMP.req1
sed 's/$/ nan/1' $TMP.req1 > $TMP.cds1
$THIS/prots_pair2dissim  -log $LOG  -power $POWER  $BLOSUM62  $AVERAGE_MODEL $TMP.req1 $TMP.univ1

cat $TMP.cds0 $TMP.cds1 $TMP.prt0 $TMP.prt1 $TMP.univ0 $TMP.univ1 | tr '\t' ' ' | sed 's|'$GENOME_DIR'/[^/]\+/[^/]\+/||g' | sed 's/\.hash-... / /g' | sed 's/\.prot-univ / /g' > $TMP.combo
$THIS/combine_dissims $TMP.combo $DISSIM_SCALE  -log $LOG > $OUT


rm -f $LOG  
rm $TMP*
