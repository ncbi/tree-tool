#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  exit 1
fi
OBJ_LIST=$1
IN_TREE=$2  # 0/1/null


if [ ! -s $OBJ_LIST ]; then
  exit 0
fi

INC=`dirname $0`


if [ $IN_TREE == 1 ]; then
  CPP_DIR/trav $OBJ_LIST "CPP_DIR/genetics/dna_mut_invert.sh $INC/../mut.dna/%h/%f $INC/../mut.index" -threads 15  -step 1 
else
  cat $OBJ_LIST >> $INC/../deleted.all
fi  


SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
BULK_REMOTE=`cat $INC/bulk_remote`

CPP_DIR/bulk.sh $SERVER $INC/bulk $BULK_REMOTE $OBJ_LIST $DATABASE..ListC

sqsh-ms  -S $SERVER  -D $DATABASE << EOT 
  update Virus
    set in_tree = $IN_TREE
    from      ListC
         join Virus on Virus.id = ListC.id
  print @@rowcount
  go -m bcp
EOT
