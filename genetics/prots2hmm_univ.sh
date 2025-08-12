#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Input: #1.prot or #1.prot_genbank"
  echo "Output: #1.{univ,prot-univ.HMM,prot-univ}"
  echo "#1: assembly file prefix"
  echo "#2: HMM library"
  echo "#3: use cut_ga (0/1)"
  echo "#4: number of cores"
  echo "#5: log file"
  exit 1
fi
PREFIX=$1
HMM_LIB=$2
CUTOFF=$3
CORES=$4
LOG=$5


IN=$PREFIX.prot_genbank
if [ ! -e $IN ]; then
  IN=$PREFIX.prot
fi


ANNOT=$PREFIX.univ
PROT_CUT=$PREFIX.prot-univ.HMM


if [ -s $IN ]; then
  TMP=$( mktemp )
  #comment $TMP

  CUTOFF_PAR=""
  if [ $CUTOFF == 1 ]; then
    CUTOFF_PAR="--cut_ga"
  fi
  hmmsearch  --tblout $TMP.hmmsearch  --domtblout $TMP.dom  --noali  -Z 10000  $CUTOFF_PAR  --cpu $CORES  $HMM_LIB $IN &>> $LOG 

  $THIS/hmmsearch2besthits $TMP.hmmsearch  -domtblout $TMP.dom  -log $LOG  > $ANNOT 
  cut -f 1,2,4,5 $ANNOT > $TMP.univ

  $THIS/filterFasta $IN  -aa  -target $TMP.univ  -replace  -cut  -len_min 20  -complexity_min 3  -log $LOG  > $PROT_CUT

  if false; then
    mkdir $TMP.seq
    $THIS/splitFasta -aa $PROT_CUT 25 $TMP.seq  -log $LOG
    $THIS/../trav $TMP.seq  -log $LOG  "hmmalign --amino --informat FASTA --outformat A2M $HMM_DIR/%f.HMM %d/%f" | sed '/^[^>]/ s/[a-z]//g' | sed '/^[[:space:]]*$/d' > $PREFIX.hmm-align
  fi

  rm -fr $TMP*
else
  cp /dev/null $ANNOT
  cp /dev/null $PROT_CUT
fi


rm -f $PREFIX.prot-univ
ln -s $( realpath $PROT_CUT ) $PREFIX.prot-univ


rm -f $LOG


