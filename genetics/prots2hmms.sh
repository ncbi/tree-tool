#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Create HMMs from proteins"
  echo "Output: files {#3,#4,#5}/(basename #1)"
  echo "#1: list of protein accessions in #2 (for seed)"
  echo "#2: protein multi-FASTA file to search against"
  echo "#3: output hmm directory"
  echo "#4: output seed directory"
  echo "#5: output stat directory"
  exit 1
fi
SEED_LIST=$1
TARGET=$2
HMM_DIR=$3
SEED_DIR=$4
STAT_DIR=$5


TMP=$( mktemp )
comment $TMP


$THIS/filterFasta $TARGET  -target $SEED_LIST  -aa > $TMP.fa

mkdir $TMP.out
NAME=$( basename $SEED_LIST )
set +o errexit
$THIS/prots2hmm.sh $TMP $TMP.fa $TARGET $TMP.out/$NAME
S=$?
set -o errexit

case $S in
  0) 
    mv $TMP.out/$NAME.stat $STAT_DIR/
    mv $TMP.out/$NAME.hmm  $HMM_DIR/
    mv $TMP.out/$NAME.SEED $SEED_DIR/
    ;;
  2) warning "Bad HMM" 
    ;;
  *) error "prots2hmm.sh exit=$S" 
    ;;
esac


rm -r $TMP*
