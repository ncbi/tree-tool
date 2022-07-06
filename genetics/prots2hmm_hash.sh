#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "$0"
  echo "#1: Protein FASTA"
  echo "#2: HMM library"
  echo "#3: Use cut_ga (0/1)"
  echo "#4: Output HMM signatitures of proteins"
  echo "#5: Output HMM hash file"
  echo "#6: log file"
  exit 1
fi
PROT=$1
HMM=$2
CUTOFF=$3
SIG=$4
HASH=$5
LOG=$6



if [ -s $PROT ]; then
  TMP=`mktemp`
  #echo $TMP

  CUT_GA=""
  if [ $CUTOFF == 1 ]; then
    CUT_GA="--cut_ga"
  fi
  hmmsearch  --tblout $TMP.hmmsearch  --noali  -Z 10000  $CUT_GA  --cpu 4  $HMM $PROT &>> $LOG
  $THIS/prots2hmm_signature $TMP.hmmsearch  -log $LOG > $SIG
  cut -f 2 $SIG | $THIS/../str2hash -log $LOG > $HASH

  rm -f $TMP*
else
  cp /dev/null $SIG
  cp /dev/null $HASH
fi

  
rm -f $LOG


