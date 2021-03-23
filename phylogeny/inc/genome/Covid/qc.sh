#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality control of distTree_inc_new.sh"
  echo "#1: go"
  echo "Requires: Time: O(n log^3(n))"
  exit 1
fi


INC=`dirname $0`

SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`


# Cf. LocusQC.sh


TMP=`mktemp`
#echo $TMP
#set -x


grep '^>' $INC/seq.fa | sed 's/^>//1' | sed 's/ .*$//1' | sort > $TMP.seq-fa

sort -u $TMP.seq-fa > $TMP.seq-fa-uniq
diff $TMP.seq-fa $TMP.seq-fa-uniq

CPP_DIR/phylogeny/tree2obj.sh $INC/tree > $TMP.tree
diff $TMP.seq-fa $TMP.tree

sqsh-ms -S $SERVER  -D $DATABASE  << EOT | sed 's/|$//1' > $TMP.locus
  select id
    from Virus
    where in_tree = 1;
  go -m bcp  
EOT
diff $TMP.seq-fa $TMP.locus

ls $INC/new/ > $TMP.new
sqsh-ms -S $SERVER  -D $DATABASE << EOT | sed 's/|$//1' | sort > $TMP.locus-new
  select id
    from Virus
    where     dead = 0
          and outlier is null
          and in_tree is null;
  go -m bcp  
EOT
#wc -l $TMP.genome-new 
#wc -l $TMP.new
diff $TMP.locus-new $TMP.new

#ls $INC/../seq/ > $TMP.seq
#diff $TMP.seq-fa $TMP.seq

echo ">aa" > $TMP.seq
echo "tttttttttttttttttttttttttt" >> $TMP.seq
blastn -query $TMP.seq  -db $INC/seq.fa | grep "Number of sequences in database:" | sed 's/,//1' | sed 's/^ *//1' > $TMP.blastn
N=`cat $TMP.seq-fa | wc -l`
M=(`cat $TMP.blastn`)
if [ $N -ne ${M[5]} ]; then
  error "$N != ${M[5]}"
fi


rm $TMP*


