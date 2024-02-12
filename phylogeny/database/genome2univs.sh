#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Input: hmm-univ.list"
  echo "#1: Genome.id"
  echo "#2: large (0/1)"
  exit 1
fi
G=$1
LARGE=$2


H=""
if [ $LARGE -eq 1 ]; then
  H=`$THIS/../../file2hash $G`
fi

$THIS/../../check_file.sh genome/$H/$G/$G.prot-univ 1


TMP=`mktemp`


grep '^>' genome/$H/$G/$G.prot-univ | sed 's/^>//1' | sed 's/ .*$//1' | sort -u > $TMP || true
N=`join  -1 1  -2 1  $TMP  hmm-univ.list | wc -l`
echo "$G $N"


rm $TMP*

