#!/bin/bash
source /home/brovervv/code/cpp/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find closest genomes"
  echo "#1: Genome.id"
  echo "#2: directory or ''"
  exit 1
fi
GENOME=$1
DIR="$2"


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
BULK_REMOTE=`cat $INC/bulk_remote`

if [ -z $DIR ]; then
  DIR=$INC/../genome
  if [ -e $INC/large ]; then
    H=`/home/brovervv/code/cpp/file2hash $GENOME`
    DIR=$DIR/$H
  fi
fi

/home/brovervv/code/cpp/database/Genome_hash_requestClosest.sh $SERVER $DATABASE $INC/bulk $BULK_REMOTE $GENOME $DIR PRT

