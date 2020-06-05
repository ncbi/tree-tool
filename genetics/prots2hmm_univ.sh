#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Input: #1.prot"
  echo "Output: #1.{univ,prot-univ,hmm-align}"
  echo "#1: assembly file prefix"
  echo "#2: HMM library"
 #echo "#3: HMM directory with .HMM-files (should match #2)"
  echo "#3: Use HMM cutoffs (0/1)"
  echo "#4: log file"
  exit 1
fi
PREFIX=$1
HMM_LIB=$2
#HMM_DIR=$3
CUTOFF=$3
LOG=$4


IN=$PREFIX.prot
ANNOT=$PREFIX.univ
PROT_CUT=$PREFIX.prot-univ


TMP=`mktemp`


CUTOFF_PAR=""
if [ $CUTOFF == 1 ]; then
  CUTOFF_PAR="--cut_ga"
fi
hmmsearch  --tblout $TMP.hmmsearch  --domtblout $TMP.dom  --noali  -Z 10000  $CUTOFF_PAR  --cpu 4  $HMM_LIB $IN &>> $LOG 

$THIS/hmmsearch2besthits $TMP.hmmsearch  -domtblout $TMP.dom  -log $LOG  > $ANNOT 
cut -f 1,2,4,5 $ANNOT > $TMP.univ

$THIS/extractFastaProt $IN $TMP.univ  -replace  -cut  -log $LOG  > $PROT_CUT

if false; then
  mkdir $TMP.seq
  $THIS/splitFastaProt $PROT_CUT 25 $TMP.seq  -log $LOG
  $THIS/../trav $TMP.seq  -log $LOG  "hmmalign --amino --informat FASTA --outformat A2M $HMM_DIR/%f.HMM %d/%f" | sed '/^[^>]/ s/[a-z]//g' | sed '/^[[:space:]]*$/d' > $PREFIX.hmm-align
fi


rm -fr $TMP*
rm -f $LOG


