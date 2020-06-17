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

INC=`dirname $0`

if [ $IN_TREE == 1 ]; then
  trav $1 "cat /home/brovervv/panfs/marker/SSU/seq/%f" >> $INC/seq.fa
else
  extractFastaDna $INC/seq.fa $OBJ_LIST  -remove > $INC/seq.fa1
  mv $INC/seq.fa1 $INC/seq.fa
fi

makeblastdb  -in $INC/seq.fa  -dbtype nucl    -logfile /dev/null


loadLISTC $OBJ_LIST

sqsh-ms  -S PROTEUS  -D uniColl << EOT 
  update Locus
    set in_tree = $IN_TREE
    from      LISTC
         join Locus on Locus.id = LISTC.id
    where     Locus.taxroot = 2759
          and Locus.gene = '18S';
  print @@rowcount;
  go -m bcp
EOT
