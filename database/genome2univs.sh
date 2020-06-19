#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: Genome.id"
  exit 1
fi
G=$1


H=`$THIS/../file2hash $G`

set +o errexit
N=`grep -c ">" genome/$H/$G/$G.prot-univ`
set -o errexit

echo "$G $N"

