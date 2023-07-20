#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control"
  echo "#1: verbose (0/1)"
  exit 1
fi
VERB=$1


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`


TMP=`mktemp`
if [ $VERB == 1 ]; then
  echo $TMP
  set -x
fi


grep '^>' $INC/prot-univ | sed 's/^>//1' | sed 's/-.*$//1' > $TMP.seq


CPP_DIR/phylogeny/tree2obj.sh $INC/tree > $TMP.tree

sort -u $TMP.seq | grep -vx -f $INC/deleted.list > $TMP.uniq
diff $TMP.uniq $TMP.tree

sqsh-ms -S $SERVER  -D $DATABASE  << EOT | sed 's/|$//1' | sort > $TMP.genome-tree
  select id
    from Prok
    where in_tree = 1
  go -m bcp  
EOT
#wc -l $TMP.genome-tree 
#wc -l $TMP.tree
diff $TMP.genome-tree $TMP.tree


CPP_DIR/phylogeny/distTree_inc_new_list.sh $INC > $TMP.new
CPP_DIR/setIntersect.sh $TMP.new $TMP.tree 0 > $TMP.intersect
if [ -s $TMP.intersect ]; then
  wc -l $TMP.intersect
  exit 1
fi

echo ">aa" > $TMP.seq_aux
echo "tttttttttttttttttttttttttt" >> $TMP.seq_aux
blastp  -query $TMP.seq_aux  -db $INC/prot-univ | grep "Number of sequences in database:" | sed 's/,//1' | sed 's/^ *//1' > $TMP.blastp
N=`cat $TMP.seq | wc -l`
M=(`cat $TMP.blastp`)
if [ $N -ne ${M[5]} ]; then
  echo "$N != ${M[5]}"
  exit 1
fi


sqsh-ms -S $SERVER  -D $DATABASE  << EOT | sed 's/|$//1' | sort > $TMP.genome-new
  select id
    from Prok
    where     outlier is null
          and in_tree is null
  go -m bcp  
EOT
diff $TMP.genome-new $TMP.new


rm $TMP*



