#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 7 ]; then
  echo "$0"
  echo "#1: Protein FASTA"
  echo "#2: HMM library"
  echo "#3: Use cut_ga (0/1)"
  echo "#4: Output HMM signatures of proteins"
  echo "#5: Output HMM hash file"
  echo "#6: number of cores"
  echo "#7: log file"
  exit 1
fi
PROT=$1
HMM=$2
CUTOFF=$3
SIG=$4
HASH=$5
CORES=$6
LOG=$7



if [ -s $PROT ]; then
  TMP=`mktemp`
 #comment $TMP

  CUT_GA=""
  if [ $CUTOFF == 1 ]; then
    CUT_GA="--cut_ga"
  fi
  hmmsearch  --tblout $TMP.hmmsearch  --noali  -Z 10000  $CUT_GA  --cpu $CORES  $HMM $PROT  &>> $LOG
  $THIS/prots2hmm_signature $TMP.hmmsearch  -log $LOG > $SIG
  cut -f 2 $SIG | $THIS/../str2hash -log $LOG > $HASH

  rm -f $TMP*
else
  > $SIG
  > $HASH
fi

  
rm -f $LOG


