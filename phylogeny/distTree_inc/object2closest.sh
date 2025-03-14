#!/bin/bash --noprofile
#source bash_common.sh
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


error "$0 is not implemented"

