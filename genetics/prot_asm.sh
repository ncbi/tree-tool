#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
FORMAT="sstart send qlen length nident"
#       1      2    3    4      5
if [ $# -ne 2 ]; then
  echo "tblastn with format: $FORMAT"
  echo "#1: Single protein FASTA file"
  echo "#2: DNA FASTA file"
  exit 1
fi
PROT=$1
DNA=$2


TMP=$( mktemp )
#comment $TMP
#set -x


GZ=$( echo $DNA | tr '.' '\n' | tail -1 )
if [ $GZ == "gz" ]; then
  gunzip -c $DNA > $TMP
  DNA=$TMP
fi

tblastn  -query $PROT  -subject $DNA  -db_gencode 11  -seg no  -comp_based_stats 0  -task tblastn-fast  -threshold 100  -window_size 15  -evalue 1e-10  -outfmt "6 $FORMAT"  1> $TMP.blast  2> /dev/null 
if [ -s $TMP.blast ]; then
  sort -k5,5nr -k3,3n $TMP.blast | head -1 
fi


rm $TMP*
