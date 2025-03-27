#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Set TypeGenome.in_tree"
  echo "#1: List of TypeGenome.id"
  echo "#2: TypeGenome.in_tree (1/null)"
  exit 1
fi
OBJ_LIST=$1
IN_TREE=$2


if [ ! -s $OBJ_LIST ]; then
  exit 0
fi


INC=$( dirname $0 )


if [ $IN_TREE == 1 ]; then
  CPP_DIR/trav $OBJ_LIST "sed 's/^>/>%f-/1' $INC/../genome/%h/%f/%f.prot-univ" >> $INC/prot-univ
else
  TMP=$( mktemp )
  CPP_DIR/trav $OBJ_LIST "CPP_DIR/genetics/fa2list.sh $INC/../genome/%h/%f/%f.prot-univ | sed 's/^/%f-/1'" > $TMP
  CPP_DIR/genetics/extractFastaProt $INC/prot-univ $TMP  -remove  -min_len 0  > $INC/prot-univ1
  rm $TMP
  mv $INC/prot-univ1 $INC/prot-univ
fi

makeblastdb  -in $INC/prot-univ  -dbtype prot    -logfile /dev/null


SERVER=$( cat $INC/server )
DATABASE=$( cat $INC/database )
BULK_REMOTE=$( cat $INC/bulk_remote )

CPP_DIR/bulk.sh $SERVER $INC/bulk $BULK_REMOTE $OBJ_LIST $DATABASE..List

sqsh-ms  -S $SERVER  -D $DATABASE << EOT 
  update TypeGenome
    set in_tree = $IN_TREE
    from      List
         join TypeGenome on TypeGenome.id = List.id
  print @@rowcount
  go -m bcp
EOT

