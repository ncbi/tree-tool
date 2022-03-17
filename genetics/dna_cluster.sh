#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Cluster DNA sequences, print representatives"
  echo "#1; directory with DNA sequence"
  echo "#2: min. fraction of identity"
  echo "#3: min. fraction of coverage"
  exit 1
fi
DIR=$1
IDENT=$2
COV=$3


TMP=`mktemp`
#echo $TMP > /dev/stderr 
#set -x


$THIS/../trav $DIR "cat %d/%f" > $TMP
$THIS/../trav $DIR "echo -e '%f\t%f'" > $TMP.pair
makeblastdb  -in $TMP   -dbtype nucl  -blastdb_version 4  -logfile /dev/null
$THIS/../trav -step 1 $DIR "blastn  -query %d/%f  -db $TMP  -dust no  -evalue 1e-10  -dbsize 10000000  -outfmt '6 qseqid sseqid qlen slen nident length qstart qend sstart send'  -num_threads 16" > $TMP.blastn
#                                                                                                                 1      2      3    4    5      6      7      8    9      10
awk '$5 >= $6 * '$IDENT' && ($8 - $7 + 1) >= $3 * '$COV' && ($10 > $9 && ($10 - $9 + 1) >= $4 * '$COV' || $10 < $9 && ($9 - $10 + 1) >= $4 * '$COV')' $TMP.blastn | cut -f 1,2 >> $TMP.pair
mkdir $TMP.dir
$THIS/../connectPairs $TMP.pair $TMP.out  -pairs 
cut -f 2 $TMP.out | sort -u


rm -r $TMP*
