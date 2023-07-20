#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find closest genomes"
  echo "#1: TypeGenome.id"
  echo "#2: directory or ''"
  exit 1
fi
GENOME=$1
DIR="$2"


INC=`dirname $0`

if [ -z $DIR ]; then
  DIR=$INC/../genome
  if [ -e $INC/large ]; then
    H=`CPP_DIR/file2hash $GENOME`
    DIR=$DIR/$H/$GENOME
  fi
fi

CPP_DIR/genetics/univ2closest.sh $DIR/$GENOME.prot-univ $INC/prot-univ

