#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Create an HMM from proteins and evaluate the HMM"
  echo "Output: files #4.{SEED,hmm,stat}"
  echo "Exit code = 2: bad HMM"
  echo "#1: temporary file name"  
  echo "#2: protein multi-FASTA file (for seed)"
  echo "#3: protein multi-FASTA file to search against"
  echo "#4: output file prefix"
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
cat $PREF.stat

if grep -q " good:1 " $PREF.stat; then
  NAME=$( basename $PREF )
  grep '^TC1: ' $PREF.stat | sed 's/^TC1: /'$NAME' /1' > $TMP.tc1
  $THIS/hmmAddCutoff $TMP.hmm $TMP.tc1 "TC" $PREF.hmm
else
  rm $PREF.SEED
  exit 2
fi
    

