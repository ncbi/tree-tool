#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 8 ]; then
  echo "Quality control of a locus incremental distance tree directory"
  echo "#1: incremental distance tree directory"
  echo "#2: SQL server name"
  echo "#3: Database"
  echo "#4: LOCUS table"
  echo "#5: LOCUS id"
  echo "#6: LOCUS.taxroot"
  echo "#7: LOCUS.gene"
  echo "#8: verbose (0/1)"
  exit 1
fi
INC=$1
SERVER=$2
DB=$3
LOCUS=$4
ID=$5
TAXROOT=$6
GENE=$7
VERB=$8


TMP=%$( mktemp )
if [ $VERB == 1 ]; then
  echo $TMP
  set -x
fi


grep '^>' $INC/seq.fa | sed 's/^>//1' | sed 's/ .*$//1' | sort > $TMP.seq-fa

sort -u $TMP.seq-fa > $TMP.seq-fa-uniq
diff $TMP.seq-fa $TMP.seq-fa-uniq

$THIS/../tree2obj.sh $INC/tree > $TMP.tree
diff $TMP.seq-fa $TMP.tree

$THIS/../distTree_inc_new_list.sh $INC > $TMP.new
$THIS/../../setIntersect.sh $TMP.new $TMP.tree 0 > $TMP.intersect
if [ -s $TMP.intersect ]; then
  error "$( wc -l $TMP.intersect )"
fi

echo ">aa" > $TMP.seq
echo "tttttttttttttttttttttttttt" >> $TMP.seq
blastn -query $TMP.seq  -db $INC/seq.fa | grep "Number of sequences in database:" | sed 's/,//g' | sed 's/^ *//1' > $TMP.blastn
N=$( < $TMP.seq-fa  wc -l )
M=( $( cat $TMP.blastn ) )
if [ $N -ne ${M[5]} ]; then
  error "$N != ${M[5]}"
fi


sqsh-ms -S $SERVER  -D $DB  -L bcp_rowsep="" << EOT | sort > $TMP.locus
  select $ID
    from $LOCUS
    where     taxroot = $TAXROOT
          and gene = '$GENE'
          and in_tree = 1;
  go -m bcp  
EOT
diff $TMP.seq-fa $TMP.locus

sqsh-ms -S $SERVER  -D $DB  -L bcp_rowsep="" << EOT | sort > $TMP.locus-new
  select $ID
    from $LOCUS
    where     taxroot = $TAXROOT
          and gene = '$GENE'
          and dead = 0
          and outlier is null
          and in_tree is null;
  go -m bcp  
EOT
diff $TMP.locus-new $TMP.new


rm $TMP*
