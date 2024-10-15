#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find closest genomes"
  echo "#1: RGenome.acc_ver"
  echo "#2: new object directory or ''"
  echo "#3: subset of Genome.id's (absolute pathname)"
  echo "#4: output file (absolute pathname)"
  exit 1
fi
GENOME=$1
DIR="$2"


exit 1
