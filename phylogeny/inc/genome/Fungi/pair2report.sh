#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print dissimilarity statistics: # identical proteins, # common universal proteins, AAI%"
  echo "#1: object 1 or new object path, where the object name is basename"
  echo "#2: #1 is new (0/1)"
  echo "#3: object 2"
  exit 1
fi
OBJ1=$1
NEW=$2
OBJ2=$3


TMP=`mktemp`
#set -x


INC=`dirname $0`
GENOME=$INC/../genome

PATH1=$OBJ1
if [ $NEW == 0 ]; then
  H1=`CPP_DIR/file2hash $OBJ1`
  PATH1=$GENOME/$H1/$OBJ1/$OBJ1
fi

H2=`CPP_DIR/file2hash $OBJ2`
PATH2=$GENOME/$H2/$OBJ2/$OBJ2


$INC/qc_object.sh $PATH1
$INC/qc_object.sh $PATH2

PRT=`CPP_DIR/setIntersect.sh $PATH1.hash-PRT $PATH2.hash-PRT 1 | wc -l`

mkdir $TMP.1
CPP_DIR/genetics/splitFasta $PATH1.prot-univ $TMP.1 -aa -noprogress
mkdir $TMP.2
CPP_DIR/genetics/splitFasta $PATH2.prot-univ $TMP.2 -aa -noprogress
CPP_DIR/trav $TMP.1 "if [ -e $TMP.2/%f ]; then CPP_DIR/dissim/seq2dissim %d/%f $TMP.2/%f -prot_name %f -blosum62 -match_len_min 20 | grep '^identity:' | cut -f 2; fi" -noprogress > $TMP.univ
UNIVS=`cat $TMP.univ | wc -l`
if [ -s $TMP.univ ]; then
  cat $TMP.univ | CPP_DIR/dm/count > $TMP.count
  AAI=`grep -w '^mean' $TMP.count | cut -f 2`
else
  AAI=0
fi
AAI_PERC=`echo "$AAI * 100" | bc`

echo -e "$PRT\t$UNIVS\t$AAI_PERC"


rm -r $TMP*