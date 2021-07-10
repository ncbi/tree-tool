#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: verbose (0/1)"
  echo "Requires: Time: O(n log^3(n))"
  exit 1
fi
VERB=$1


INC=`dirname $0`

SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`


# Cf. LocusQC.sh


CPP_DIR/phylogeny/distTree_inc_indiscern_qc.sh $INC
sort -cu $INC/good


TMP=`mktemp`
if [ $VERB == 1 ]; then
  echo $TMP
  set -x
fi


CPP_DIR/phylogeny/tree2obj.sh $INC/tree > $TMP.tree


sqsh-ms -S $SERVER  -D $DATABASE  << EOT | sed 's/|$//1' | sort > $TMP.db-tree
  select id
    from Virus
    where in_tree = 1
  go -m bcp  
EOT
diff $TMP.tree $TMP.db-tree

CPP_DIR/phylogeny/distTree_inc_new_list.sh $INC > $TMP.new
sqsh-ms -S $SERVER  -D $DATABASE << EOT | sed 's/|$//1' | sort > $TMP.db-new
  select id
    from Virus
    where     dead = 0
          and outlier is null
          and in_tree is null
  go -m bcp  
EOT
#wc -l $TMP.genome-new 
#wc -l $TMP.new
diff $TMP.db-new $TMP.new

CPP_DIR/setIntersect.sh $TMP.new $TMP.tree 0 > $TMP.inter
if [ -s $TMP.inter ]; then
  error "$INC/new are in the tree"
fi

sort -u $INC/../deleted.all > $TMP.deleted

CPP_DIR/setIntersect.sh $TMP.deleted $TMP.tree 0 > $TMP.inter
if [ -s $TMP.inter ]; then
  error "deleted.all are in the tree"
fi

CPP_DIR/setIntersect.sh $TMP.deleted $TMP.new 0 > $TMP.inter
if [ -s $TMP.inter ]; then
  error "deleted.all are in the $INC/new"
fi

sort -u $INC/good > $TMP.good

CPP_DIR/setIntersect.sh $TMP.deleted $TMP.good 0 > $TMP.inter
if [ -s $TMP.inter ]; then
  error "deleted.all are in the $INC/good"
fi

sqsh-ms -S $SERVER  -D $DATABASE << EOT | sed 's/|$//1' | sort > $TMP.bad
  select id
    from Virus
    where    dead = 1
          or outlier is not null
  go -m bcp  
EOT

CPP_DIR/setIntersect.sh $TMP.bad $TMP.good 0 > $TMP.inter
if [ -s $TMP.inter ]; then
  error "Bad are in the $INC/good"
fi


rm $TMP*


