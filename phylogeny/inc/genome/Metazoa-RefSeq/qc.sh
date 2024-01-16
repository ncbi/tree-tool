#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: verbose (0/1)"
  exit 1
fi
VERB=$1


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`


if [ -e $INC/good ]; then
  sort -c -u $INC/good
fi


TMP=`mktemp`
if [ $VERB == 1 ]; then
  echo $TMP
  set -x
fi


CPP_DIR/phylogeny/tree2obj.sh $INC/tree > $TMP.tree
sqsh-ms -S $SERVER  -D $DATABASE  << EOT | sed 's/|$//1' | sort > $TMP.genome-tree
  select acc_ver
    from RGenome
    where     tax_root = 33208
          and in_tree = 1
  go -m bcp  
EOT
#wc -l $TMP.genome-tree 
#wc -l $TMP.tree
diff $TMP.genome-tree $TMP.tree


CPP_DIR/phylogeny/distTree_inc_new_list.sh $INC > $TMP.new
sqsh-ms -S $SERVER  -D $DATABASE  << EOT | sed 's/|$//1' | sort > $TMP.genome-new
  select acc_ver
    from RGenome
    where     tax_root = 33208
          and dead = 0
          and prots is not null
          and outlier is null
          and in_tree is null
  go -m bcp  
EOT
diff $TMP.genome-new $TMP.new


rm $TMP*
