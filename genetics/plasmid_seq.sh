#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Create files for #1 in the current directory"
  echo "Input: ./#1"
  echo "Invokes: prodigal"
  echo "#1: plasmid accession"
 #echo "#2: Pfam HMM library"
  echo "#2: iniveral HMM library"
  echo "#3: list of ribosomal HMMs"
  echo "#4: log file"
  echo "Time: 40 sec."
  exit 1
fi
PL=$1
#PFAM=$2
UNIV=$2
RIB=$3
LOG=$4


if [ ! -e $PL ]; then
  echo "No $PL sequence" >> $LOG
  exit 1
fi
if [ ! -s $PL ]; then
  echo "Empty genome" >> $LOG
  exit 1
fi

$THIS/dna2stat $PL  -log $LOG > $PL.stat

prodigal  -o /dev/null  -i $PL  -d $PL.cds  -c  -p meta  &>> $LOG

if [ -s $PL.cds ]; then
 #gunzip $PL.prot.gz
  # $PL.prot
  $THIS/orf2prot $PL.cds  -log $LOG > $PL.prot
  # Find mobilome (transposases) ??
  
  $THIS/prots2hmm_univ.sh $PL $UNIV 1 4 $LOG
  grep -w -f $RIB $PL.univ > $PL.ribosomal || true
  
 #$THIS/prots2hmm_hash.sh $PL.prot $PFAM 0 $PL.HMM $PL.hash-HMM 4 $LOG
 #cut -f 2 $PL.HMM | tr ',' '\n' | sort -u > $PL.pfam
 #gzip $PL.HMM
  
  gzip $PL.prot
fi

rm $PL.cds



