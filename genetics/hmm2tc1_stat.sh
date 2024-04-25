#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print .tsv-file with TC1/LENG statistics for an HMM library"
  echo "#1: HMM library"
  exit 1
fi
LIB=$1


TMP=`mktemp`
#comment $TMP
#set -x


mkdir $TMP.hmm
$THIS/hmmSplit $LIB $TMP.hmm -qc

mkdir $TMP.pair
$THIS/../trav $TMP.hmm "grep -w %q^LENG\|^TC%q %d/%f | awk %q{OFS=%Q%bt%Q;print %D1, %D2};%q> $TMP.pair/%f" -step 1

$THIS/../tsv/pairs2tsv $TMP.pair  -file_col "hmm"  -qc 


rm -r $TMP*
