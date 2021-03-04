#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Test k-mer index"
  echo "#1: go"
  exit 1
fi


TMP=`mktemp`
echo $TMP


$THIS/kmerIndex_make $TMP.kmi 14 -qc
$THIS/kmerIndex_add  $TMP.kmi $THIS/data/5_8S.fa  -qc
$THIS/kmerIndex_stat $TMP.kmi -qc

$THIS/fa2list.sh $THIS/data/5_8S.fa | sort -R > $TMP.list-all
head $TMP.list-all | sort > $TMP.list
$THIS/extractFastaDna $THIS/data/5_8S.fa $TMP.list -qc > $TMP.fa
mkdir $TMP.seq
$THIS/splitFasta $TMP.fa $TMP.seq -qc 

mkdir $TMP.out
$THIS/../trav $TMP.seq "$THIS/kmerIndex_find $TMP.kmi %d/%f 100 -qc -self > $TMP.out/%f"
$THIS/../trav $TMP.out '[ `cat %d/%f | wc -l` == 100 ]'
$THIS/../trav $TMP.out 'grep -w %f %d/%f'


rm -r $TMP*
