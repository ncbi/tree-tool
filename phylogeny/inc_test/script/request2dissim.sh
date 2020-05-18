#!/bin/bash
THIS=`dirname $0`
if [ $# -ne 3 ]; then
  echo "Compute dissimilarities"
  echo "#1: file with pairs <Object1> <Object2>"
  echo "#2: output file"
  echo "#3: log (temporary)"
  exit 1
fi
REQ=$1
OUT=$2
LOG=$3


TMP=`mktemp`
#echo $TMP 


HASH_DIR=$THIS/../hash


cat $REQ | awk '{print $1};' > $TMP.1
cat $REQ | awk '{print $2};' > $TMP.2
DT/trav $TMP.1 "echo %n $HASH_DIR/%f/%f" > $TMP.f1
DT/trav $TMP.2 "echo %n $HASH_DIR/%f/%f" > $TMP.f2
join  -1 1  -2 1  $TMP.f1 $TMP.f2 | cut -d ' ' -f 2,3 > $TMP.req


DT/phylogeny/hash_request2dissim $TMP.req  $OUT  -intersection_min 50  -ratio_min 0.5  -log $LOG 


rm $TMP*
