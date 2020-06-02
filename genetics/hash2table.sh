#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "#1: uniColl..Genome.id"
  echo "#2: hash type (uniColl..GenomeHash.type)"
  exit 1
fi
GENOME=$1
TYPE=$2


H=`$THIS/../file2hash $GENOME`
cat genome/$H/$GENOME/$GENOME.hash-$TYPE | sed 's/^/'$GENOME'#'$TYPE'#/1' | tr '#' '\t'


