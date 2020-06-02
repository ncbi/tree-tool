#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Create an HMM from proteins and evaluate the HMM"
  echo "Output: files #4.{SEED,hmm,stat}"
  echo "Exit code = 2: bad HMM"
  echo "#1: Temporary file name"  
  echo "#2: Protein file in FASTA format (for seed)"
  echo "#3: Protein file in FASTA format to search against"
  echo "#4: Output file prefix"
  exit 1
fi
TMP=$1
SEED=$2
TARGET=$3
PREF=$4


if [ ! -e $PREF.SEED ]; then
  # PAR
 #clustalw -infile=$2 -outfile=$4.SEED -output=fasta -quicktree > /dev/null
  muscle -sv -noanchors  -distance1 kbit20_3  -maxiters 2  -maxtrees 99  -in $SEED  -out $PREF.SEED >& /dev/null
    # >1 alignment ??
fi

# PAR
hmmbuild  --informat afa  $TMP.hmm $PREF.SEED > /dev/null  

# No TC
hmmsearch  --tblout $TMP.hmmsearch  --noali  -Z 10000  $TMP.hmm $TARGET >& /dev/null 

# PAR
# GA: must return neighboring families ??
$THIS/hmm_tc1 -noprogress  -hmmsearch $TMP.hmmsearch  -error_max 0.01  -seqFName $SEED > $PREF.stat
cat $OREF.stat

set -o errexit
grep " good:1 " $PREF.stat > /dev/null
S=$?
set -o errexit
if [ $S == 0 ]; then
  grep '^TC1: ' $PREF.stat > $TMP.tc1
  $THIS/hmmAddTC1 $TMP.hmm $TMP.tc1 $PREF.hmm
else
  rm $PREF.SEED
  exit 2
fi
    

