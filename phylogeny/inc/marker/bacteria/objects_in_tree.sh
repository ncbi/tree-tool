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

BASE_DIR=`dirname $0`

if [ $IN_TREE == 1 ]; then
  CPP_DIR/trav $1 "cat $BASE_DIR/../seq/%f" >> $BASE_DIR/seq.fa
else
  CPP_DIR/genetics/extractFastaDna $BASE_DIR/seq.fa $OBJ_LIST  -remove > $BASE_DIR/seq.fa1
  mv $BASE_DIR/seq.fa1 $BASE_DIR/seq.fa
fi

makeblastdb  -in $BASE_DIR/seq.fa  -dbtype nucl  -blastdb_version 4  -logfile /dev/null

if false; then
loadLISTC $OBJ_LIST
sqsh-ms  -S ""  -D uniColl << EOT 
  update Locus
    set in_tree = $IN_TREE
    from      LISTC
         join Locus on Locus.id = LISTC.id;
  print @@rowcount;
  go -m bcp
EOT
fi
