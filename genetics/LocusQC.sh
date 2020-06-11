#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Quality control of a locus incremental distance tree directory"
  echo "#1: incremental distance tree directory"
  echo "#2: SQL server name or ''"
  echo "#3: Database"
  echo "#4: LOCUS table"
  echo "#5: LOCUS.taxroot"
  echo "#6: LOCUS.gene"
  exit 1
fi
INC=$1
SERVER="$2"
DB=$3
LOCUS=$4
TAXROOT=$5
GENE=$6


TMP=`mktemp`
#echo $TMP
#set -x


grep '^>' $INC/seq.fa | sed 's/^>//1' | sed 's/ .*$//1' | sort > $TMP.seq-fa

sort -u $TMP.seq-fa > $TMP.seq-fa-uniq
diff $TMP.seq-fa $TMP.seq-fa-uniq

$THIS/../phylogeny/tree2obj.sh $INC/tree > $TMP.tree
diff $TMP.seq-fa $TMP.tree

ls $INC/new/ > $TMP.new
$THIS/../setIntersect.sh $TMP.new $TMP.tree 0 > $TMP.intersect
if [ -s $TMP.intersect ]; then
  wc -l $TMP.intersect
  exit 1
fi

echo ">aa" > $TMP.seq
echo "tttttttttttttttttttttttttt" >> $TMP.seq
blastn -query $TMP.seq  -db $INC/seq.fa | grep "Number of sequences in database:" | sed 's/,//1' | sed 's/^ *//1' > $TMP.blastn
N=`cat $TMP.seq-fa | wc -l`
M=(`cat $TMP.blastn`)
if [ $N -ne ${M[5]} ]; then
  echo "$N != ${M[5]}"
  exit 1
fi


if [ $SERVER ]; then
sqsh-ms -S $SERVER  -D $DB  << EOT | sed 's/|$//1' > $TMP.locus
  select id
    from $LOCUS
    where     taxroot = $TAXROOT
          and gene = '$GENE'
          and in_tree = 1;
  go -m bcp  
EOT
diff $TMP.seq-fa $TMP.locus

sqsh-ms -S $SERVER  -D $DB   << EOT | sed 's/|$//1' | sort > $TMP.locus-new
  select id
    from $LOCUS
    where     taxroot = $TAXROOT
          and gene = '$GENE'
          and dead = 0
          and outlier is null
          and in_tree is null;
  go -m bcp  
EOT
#wc -l $TMP.genome-new 
#wc -l $TMP.new
diff $TMP.locus-new $TMP.new
fi


rm $TMP*
