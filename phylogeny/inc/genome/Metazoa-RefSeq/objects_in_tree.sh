#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Set RGenome.in_tree"
  echo "#1: List of RGenome.acc_ver"
  echo "#2: RGenome.in_tree (1/null)"
  exit 1
fi
OBJ_LIST=$1
IN_TREE=$2


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
BULK_REMOTE=`cat $INC/bulk_remote`

CPP_DIR/bulk.sh $SERVER $INC/bulk $BULK_REMOTE $OBJ_LIST $DATABASE..ListC

sqsh-ms  -S $SERVER  -D $DATABASE << EOT 
  update RGenome
    set in_tree = $IN_TREE
    from      ListC
         join RGenome on RGenome.acc_ver = ListC.id
  print @@rowcount
  go -m bcp
EOT

