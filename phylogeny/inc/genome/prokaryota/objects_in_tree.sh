#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Set Prok.in_tree"
  echo "#1: List of Prok.id"
  echo "#2: Prok.in_tree (1/null)"
  exit 1
fi
OBJ_LIST=$1
IN_TREE=$2


INC=`dirname $0`


if [ $IN_TREE == 1 ]; then
  CPP_DIR/trav $OBJ_LIST "sed 's/^>/>%f-/1' $INC/../genome/%f/%f.prot-univ" >> $INC/prot-univ
  makeblastdb  -in $INC/prot-univ  -dbtype prot    -logfile /dev/null
else
  cat $OBJ_LIST >> $INC/deleted.list
fi


SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
BULK_REMOTE=`cat $INC/bulk_remote`

CPP_DIR/bulk.sh $SERVER $INC/bulk $BULK_REMOTE $OBJ_LIST $DATABASE..List

sqsh-ms  -S $SERVER  -D $DATABASE << EOT 
  update Prok
    set in_tree = $IN_TREE
    from      List
         join Prok on Prok.id = List.id
  print @@rowcount
  go -m bcp
EOT

