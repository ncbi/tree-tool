#!/bin/bash 
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
  CPP_DIR/trav $1 "cat $INC/../seq/%f" >> $INC/seq.fa
else
  CPP_DIR/genetics/extractFastaDna $INC/seq.fa $OBJ_LIST  -remove > $INC/seq.fa1
  mv $INC/seq.fa1 $INC/seq.fa
fi

makeblastdb  -in $INC/seq.fa  -dbtype nucl    -logfile /dev/null


SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
BULK_REMOTE=`cat $INC/bulk_remote`

CPP_DIR/database/bulk.sh $SERVER $INC/bulk $BULK_REMOTE $OBJ_LIST $DATABASE..ListC

sqsh-ms  -S $SERVER  -D $DATABASE << EOT 
  update Virus
    set in_tree = $IN_TREE
    from      ListC
         join Virus on Virus.id = ListC.id;
  print @@rowcount;
  go -m bcp
EOT