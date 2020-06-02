#!/bin/bash 
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Output: genome/#1/, hash-PRT/#1, hash-CDS/#1"
  echo "Input: genome/#1/#1.dna: delete"
  echo "#1: PDT.id"
  echo "#2: MLST scheme taxid | 0"
  echo "#3: log file (absolute path)"
  exit 1
fi
PDT=$1
MLST=$2
LOG=$3


cp /dev/null $LOG


cd genome/$PDT/

if [ ! -e $PDT.dna ]; then
  echo "No $PDT.dna" >> $LOG
  exit 1
fi
if [ ! -s $PDT.dna ]; then
  echo "Empty genome" >> $LOG
  exit 1
fi

$THIS/prok_seq.sh $PDT $PDT.dna null $MLST $LOG 0

cd ../.. 


DIR=`pwd`

rm -f hash-CDS/$PDT
ln -s $DIR/genome/$PDT/$PDT.hash-CDS hash-CDS/$PDT

rm -f hash-PRT/$PDT
ln -s $DIR/genome/$PDT/$PDT.hash-PRT hash-PRT/$PDT


rm genome/$PDT/$PDT.dna


rm -f $LOG


