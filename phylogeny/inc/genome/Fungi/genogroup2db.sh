#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "$0"
  echo "#1: File genogroup_table"
  exit 1
fi
TAB=$1


INC=$( dirname $0 )
CPP_DIR/phylogeny/database/genogroup2db.sh $INC $TAB
