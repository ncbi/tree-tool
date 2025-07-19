#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Quality control of #1"
  echo "#1: inc/"
  echo "#2: verbose (0/1)"
  exit 1
fi
INC=$1
VERB=$2


SERVER=$( cat $INC/server )
DATABASE=$( cat $INC/database )
TAXROOT=$( cat $INC/../taxroot )


if [ -e $INC/good ]; then
  sort -c -u $INC/good
fi


TMP=$( mktemp )
if [ $VERB == 1 ]; then
  echo $TMP
  set -x
fi


$THIS/../tree2obj.sh $INC/tree > $TMP.tree
sqsh-ms -S $SERVER  -D $DATABASE  -L bcp_rowsep="" << EOT | sort > $TMP.genome-tree
  select genome
    from GenomeTaxroot
    where     taxroot = $TAXROOT
          and in_tree = 1
  go -m bcp  
EOT
#wc -l $TMP.genome-tree 
#wc -l $TMP.tree
diff $TMP.genome-tree $TMP.tree


$THIS/../distTree_inc_new_list.sh $INC > $TMP.new
sqsh-ms -S $SERVER  -D $DATABASE  -L bcp_rowsep=""  << EOT | sort > $TMP.genome-new
  select id
    from GenomeTaxroot_vw
    where     taxroot = $TAXROOT
          and dead = 0
          and prots is not null
          and hybrid is null
          and outlier is null
          and in_tree is null
  go -m bcp  
EOT
diff $TMP.genome-new $TMP.new


rm $TMP*
