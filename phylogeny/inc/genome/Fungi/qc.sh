#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: verbose (0/1)"
  exit 1
fi
VERB=$1


INC=$( dirname $0 )
CPP_DIR/phylogeny/database/GenomeQC.sh $INC $VERB

