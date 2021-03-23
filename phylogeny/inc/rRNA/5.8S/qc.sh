#!/bin/bash
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: go"
  echo "Requires: Time: O(n log^3(n))"
  exit 1
fi


INC=`dirname $0`
SERVER=`cat $INC/server`
DB=`cat $INC/database`


TMP=`mktemp`
#echo $TMP
#set -x


CPP_DIR/genetics/kmerIndex_stat $INC/seq.kmi -qc > $TMP.stat
N_KMI=`grep '# DNA sequences:' $TMP.stat | tail -1 | sed 's/^.*: //1'`

CPP_DIR/phylogeny/tree2obj.sh $INC/tree > $TMP.tree
N_TREE=`cat $TMP.tree | wc -l`
if [ $N_KMI != $N_TREE ]; then
  error "K-mer index size ($N_KMI) != tree size ($N_TREE)"
fi

CPP_DIR/phylogeny/distTree_inc_new_list.sh $INC > $TMP.new
CPP_DIR/setIntersect.sh $TMP.new $TMP.tree 0 > $TMP.intersect
if [ -s $TMP.intersect ]; then
  wc -l $TMP.intersect
  exit 1
fi


sqsh-ms -S $SERVER  -D $DB  << EOT | sed 's/|$//1' | sort > $TMP.locus
  select accession
    from Locus5_8S
    where     selected = 1
          and in_tree = 1;
  go -m bcp  
EOT
diff $TMP.tree $TMP.locus

sqsh-ms -S $SERVER  -D $DB  << EOT | sed 's/|$//1' | sort > $TMP.locus-new
  select accession
    from Locus5_8S
    where     selected = 1
          and in_tree is null
          and outlier is null;
  go -m bcp  
EOT
diff $TMP.locus-new $TMP.new


rm $TMP*
