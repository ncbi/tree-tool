#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "$0"
  echo "#1: File genogroup_table"
  exit 1
fi
IN=$1


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
BULK_REMOTE=`cat $INC/bulk_remote`

CPP_DIR/bulk.sh $SERVER $INC/bulk $BULK_REMOTE $IN $DATABASE..ListC


sqsh-ms  -S $SERVER  -D $DATABASE  -L exit_failcount=1 << EOF | sed 's/|$//1' 
  EXEC Genogroup2outliers 2, 1
  go -m bcp
EOF

