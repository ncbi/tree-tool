#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Test k-mer index"
  echo "#1: go"
  exit 1
fi


TMP=$( mktemp )
comment $TMP


section "make"
$THIS/kmerIndex_make $TMP.kmi 14 -qc

section "add"
$THIS/kmerIndex_add  $TMP.kmi $THIS/data/5_8S.fa  -qc

section "stat"
$THIS/kmerIndex_stat $TMP.kmi -qc

section "$TMP.seq/"
$THIS/fa2list.sh $THIS/data/5_8S.fa > $TMP.acc
$THIS/../setRandOrd $TMP.acc -sigpipe -qc | head -10 | sort > $TMP.list
#$THIS/extractFastaDna $THIS/data/5_8S.fa $TMP.list -qc > $TMP.fa
$THIS/filterFasta $THIS/data/5_8S.fa  -target $TMP.list  -qc > $TMP.fa
mkdir $TMP.seq
$THIS/splitFasta $TMP.fa $TMP.seq -qc 

section "find"
mkdir $TMP.out
$THIS/../trav $TMP.seq "$THIS/kmerIndex_find $TMP.kmi %d/%f 100 -qc -self > $TMP.out/%f"
$THIS/../trav $TMP.out 'grep -w %f %d/%f'
$THIS/../trav $TMP.out 'cat %d/%f | wc -l | grep -w 100'


rm -r $TMP*
