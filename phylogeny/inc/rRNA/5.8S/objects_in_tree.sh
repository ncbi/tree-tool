#!/bin/bash  --noprofile
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
KMER=14

if [ $IN_TREE == 1 ]; then
  CPP_DIR/trav $OBJ_LIST "cat $INC/../seq/%h/%f" > $INC/seq.fa
  if [ ! -e $INC/seq.kmi ]; then
    CPP_DIR/genetics/kmerIndex_make $INC/seq.kmi $KMER -qc
  fi
  CPP_DIR/genetics/kmerIndex_add  $INC/seq.kmi $INC/seq.fa -qc
else
  CPP_DIR/phylogeny/tree2obj.sh $INC/tree > $INC/tree.list
  CPP_DIR/setMinus $INC/tree.list $OBJ_LIST > $INC/tree.list-new
  rm $INC/tree.list
  wc -l $INC/tree.list-new
  CPP_DIR/trav $INC/tree.list-new "cat $INC/../seq/%h/%f" > $INC/seq.fa
  rm $INC/tree.list-new
  CPP_DIR/genetics/kmerIndex_make $INC/seq.kmi $KMER -qc
  CPP_DIR/genetics/kmerIndex_add  $INC/seq.kmi $INC/seq.fa -qc
fi
rm $INC/seq.fa


SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
BULK_REMOTE=`cat $INC/bulk_remote`

CPP_DIR/bulk.sh $SERVER $INC/bulk $BULK_REMOTE $OBJ_LIST $DATABASE..ListC

sqsh-ms  -S $SERVER  -D $DATABASE << EOT 
  update Locus5_8S
    set in_tree = $IN_TREE
    from      ListC
         join Locus5_8S on     Locus5_8S.accession = ListC.id
                           and Locus5_8S.selected = 1
  print @@rowcount
  go -m bcp
EOT
