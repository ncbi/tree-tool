#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print: <Genome.id> <univs>"
  echo "Input: hmm-univ.list"
  echo "#1: Genome.id"
  echo "#2: large (0/1)"
  echo "#3: suffix: HMM|tblastn"
  exit 1
fi
G=$1
LARGE=$2
SUF=$3


sort -cu hmm-univ.list


H=""
if [ $LARGE -eq 1 ]; then
  H=$( $THIS/../../file2hash $G )
fi

F=genome/$H/$G/$G.prot-univ.$SUF
if [ ! -e $F ]; then
  exit 0
fi


TMP=$( mktemp )


grep '^>' $F | sed 's/^>//1' | sed 's/ .*$//1' | sort -u > $TMP || true
N=$( join  -1 1  -2 1  $TMP  hmm-univ.list | wc -l )
echo "$G $N"


rm $TMP*

