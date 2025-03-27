#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control"
  echo "#1: verbose (0/1)"
  exit 1
fi
VERB=$1


INC=$( dirname $0 )
SERVER=$( cat $INC/server )
DATABASE=$( cat $INC/database )


TMP=$( mktemp )
if [ $VERB == 1 ]; then
  echo $TMP
  set -x
fi


CPP_DIR/phylogeny/tree2obj.sh $INC/tree > $TMP.tree
sqsh-ms -S $SERVER  -D $DATABASE  << EOT | sed 's/|$//1' | sort > $TMP.genome-tree
  select id
    from TypeGenome
    where in_tree = 1
  go -m bcp  
EOT
#wc -l $TMP.genome-tree 
#wc -l $TMP.tree
diff $TMP.genome-tree $TMP.tree


CPP_DIR/phylogeny/distTree_inc_new_list.sh $INC > $TMP.new
sqsh-ms -S $SERVER  -D $DATABASE  << EOT | sed 's/|$//1' | sort > $TMP.genome-new
  select id
    from TypeGenome
    where     outlier is null
          and in_tree is null
  go -m bcp  
EOT
diff $TMP.genome-new $TMP.new


CPP_DIR/genetics/fa2list.sh $INC/prot-univ | sort > $TMP.seq
sort -u $TMP.seq > $TMP.seq-uniq
diff $TMP.seq $TMP.seq-uniq

sed 's/-.*$//1' $TMP.seq | sort -u > $TMP.seq-list
diff $TMP.seq-list $TMP.tree


echo ">aa" > $TMP.aa
echo "GKRVLVMGLGLQGSGMAAARYAAQQGAIVRVTDMKSPDILAPSVRALAGLPIEFILGQHREEDFIWADIVIRNPGVPRTS" >> $TMP.aa
blastp -query $TMP.aa  -db $INC/prot-univ | grep "Number of sequences in database:" | sed 's/,//g' | sed 's/^ *//1' > $TMP.blastp
N=$( < $TMP.seq  wc -l )
M=( $( cat $TMP.blastp ) )
if [ $N -ne ${M[5]} ]; then
  error "$N != ${M[5]}"
fi


rm $TMP*
