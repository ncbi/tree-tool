#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find closest objects"
  echo "#1: object"
  echo "#2: directory or ''"
  echo "#3: subset of Genome.id's (absolute pathname)"
  echo "#4: output file (absolute pathname)"
  exit 1
fi
OBJ=$1
DIR="$2"


error "$0 is not implemented"

