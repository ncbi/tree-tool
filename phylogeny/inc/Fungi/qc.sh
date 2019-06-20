#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: Incremental distance tree directory"
  echo "Requires: Time: O(n log^3(n))"
  exit 1
fi
INC=$1


TMP=`mktemp`
#echo $TMP


tree2obj.sh $INC/tree > $TMP.tree
sqsh-ms -S PROTEUS  -D uniColl  -U anyone  -P allowed  << EOT | sed 's/|$//1' | sort > $TMP.genome-tree
  select id
    from Genome
    where     tax_root = 4751
          and in_tree = 1;
  go -m bcp  
EOT
diff $TMP.genome-tree $TMP.tree


ls $INC/new/ > $TMP.new
sqsh-ms -S PROTEUS  -D uniColl  -U anyone  -P allowed  << EOT | sed 's/|$//1' | sort > $TMP.genome-new
  select id
    from Genome
    where     tax_root = 4751
          and dead = 0
          and prots is not null
          and outlier is null
          and in_tree is null;
  go -m bcp  
EOT
diff $TMP.genome-new $TMP.new


rm $TMP*
