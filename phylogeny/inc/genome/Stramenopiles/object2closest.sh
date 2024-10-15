#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find closest genomes among #3"
  echo "#1: Genome.id"
  echo "#2: new object directory or ''"
  echo "#3: subset of Genome.id's (absolute pathname)"
  echo "#4: output file (absolute pathname)"
  exit 1
fi
GENOME=$1
DIR="$2"
SUBSET=$3
OUT=$4


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
BULK_REMOTE=`cat $INC/bulk_remote`
TAX=`cat $INC/../tax_id`

if [ -z $DIR ]; then
  DIR=$INC/../genome
  if [ -e $INC/large ]; then
    H=`CPP_DIR/file2hash $GENOME`
    DIR=$DIR/$H/$GENOME
  fi
fi

CPP_DIR/phylogeny/database/Genome_hash_requestClosest.sh $SERVER $DATABASE $INC/bulk $BULK_REMOTE $GENOME $TAX $DIR "PRT" $SUBSET > $OUT



