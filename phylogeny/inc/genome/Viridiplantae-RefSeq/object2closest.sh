#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find closest genomes"
  echo "#1: RGenome.acc_ver"
  echo "#2: new object directory or ''"
  exit 1
fi
GENOME=$1
DIR="$2"


exit 1
