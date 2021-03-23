#!/bin/bash 
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  exit 1
fi
OBJ_LIST=$1
IN_TREE=$2  # 0/1/null


if [ $IN_TREE == 0 ]; then
  error "Cannot remove sequences from a K-mer index"
fi

if [ ! -s $OBJ_LIST ]; then
  exit 0
fi

INC=`dirname $0`

CPP_DIR/trav $OBJ_LIST "cat $INC/../seq/%h/%f" > $INC/seq.fa
if [ ! -e $INC/seq.kmi ]; then
  CPP_DIR/genetics/kmerIndex_make $INC/seq.kmi 14 -qc
fi
CPP_DIR/genetics/kmerIndex_add  $INC/seq.kmi $INC/seq.fa -qc
rm $INC/seq.fa


SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
BULK_REMOTE=`cat $INC/bulk_remote`

CPP_DIR/database/bulk.sh $SERVER $INC/bulk $BULK_REMOTE $OBJ_LIST $DATABASE..ListC

sqsh-ms  -S $SERVER  -D $DATABASE << EOT 
  update Locus5_8S
    set in_tree = $IN_TREE
    from      ListC
         join Locus5_8S on     Locus5_8S.accession = ListC.id
                           and Locus5_8S.selected = 1;
  print @@rowcount;
  go -m bcp
EOT
