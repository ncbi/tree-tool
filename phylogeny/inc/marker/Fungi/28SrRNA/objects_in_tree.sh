#!/bin/bash 
source bash_common.sh
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
  trav $1 "cat /home/brovervv/panfs/marker/Fungi/28SrRNA/seq/%f" >> $BASE_DIR/seq.fa
else
  extractFastaDna $BASE_DIR/seq.fa $OBJ_LIST  -remove > $BASE_DIR/seq.fa1
  mv $BASE_DIR/seq.fa1 $BASE_DIR/seq.fa
fi

makeblastdb  -in $BASE_DIR/seq.fa  -dbtype nucl  -blastdb_version 4  -logfile /dev/null


loadLISTC $OBJ_LIST

sqsh-ms  -S PROTEUS  -D uniColl << EOT 
  update Locus
    set in_tree = $IN_TREE
    from      LISTC
         join Locus on Locus.id = LISTC.id;
  print @@rowcount;
  go -m bcp
EOT
