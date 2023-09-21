#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Input: hmm-univ.list"
  echo "#1: Genome.id"
  exit 1
fi
G=$1


H=`$THIS/../../file2hash $G`

if false; then
  set +o errexit
  N=`grep -c ">" genome/$H/$G/$G.prot-univ`
  set -o errexit
fi


TMP=`mktemp`


grep '^>' genome/$H/$G/$G.prot-univ | sed 's/^>//1' | sed 's/ .*$//1' | sort -u > $TMP || true
N=`join  -1 1  -2 1  $TMP  hmm-univ.list | wc -l`
echo "$G $N"


rm $TMP*

