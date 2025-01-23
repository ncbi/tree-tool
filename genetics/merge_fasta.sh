#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Append sequences in #1 to FASTA files named as the sequence identifiers"
  echo "#1: FASTA"
  echo "#2: aa|na"
  echo "#3: name to be appended to sequence identifiers"
  echo "#4: output directory"
  exit 1
fi
IN=$1
AA=$2
NAME=$3
DIR=$4


if [ -z "$NAME" ]; then
  error "Empty name"
fi


TMP=$( mktemp )
#comment $TMP 
#set -x


mkdir $TMP.dir
$THIS/splitFasta $IN $TMP.dir  -$AA -pseudo -whole  -suffix $NAME  -noprogress
$THIS/../trav $TMP.dir "cat %d/%f >> $DIR/%f"  -noprogress


rm -r $TMP*
