#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Cluster DNA sequences, if #1 = #2 then print representatives else print #1-#2 pairs"
  echo "#1; query directory with DNA sequence"
  echo "#2; subject directory with DNA sequence"
  echo "#3: min. fraction of identity"
  echo "#4: min. fraction of coverage"
  exit 1
fi
QDIR=$1
SDIR=$2
IDENT=$3
COV=$4


TMP=`mktemp`
#echo $TMP > /dev/stderr 
#set -x


$THIS/../trav $SDIR "cat %d/%f" > $TMP
makeblastdb  -in $TMP   -dbtype nucl  -blastdb_version 4  -logfile /dev/null

# $TMP.pair
if [ $QDIR == $SDIR ]; then
  $THIS/../trav $QDIR "echo -e '%f\t%f'" > $TMP.pair
fi
$THIS/../trav -step 1 $QDIR "blastn  -query %d/%f  -db $TMP  -dust no  -evalue 1e-10  -dbsize 10000000  -max_target_seqs 100000 -outfmt '6 qseqid sseqid qlen slen nident length qstart qend sstart send'  -num_threads 16" > $TMP.blastn
#                                                                                                                                          1      2      3    4    5      6      7      8    9      10
awk '$5 >= $6 * '$IDENT' && ($8 - $7 + 1) >= $3 * '$COV' && ($10 > $9 && ($10 - $9 + 1) >= $4 * '$COV' || $10 < $9 && ($9 - $10 + 1) >= $4 * '$COV')' $TMP.blastn | cut -f 1,2 >> $TMP.pair

if [ $QDIR == $SDIR ]; then
  $THIS/../connectPairs $TMP.pair $TMP.out  -pairs 
  cut -f 2 $TMP.out | sort -u
else
  cat $TMP.pair
fi


rm $TMP*
