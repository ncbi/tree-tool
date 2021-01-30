#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Symmetric best hits dissimilarity by BLASTP"
  echo "#1: gzip'ed protein FASTA 1"
  echo "#2: gzip'ed protein FASTA 2"
  exit 1
fi
PROTGZ1=$1
PROTGZ2=$2


TMP=`mktemp`
#echo $TMP 


gunzip $PROTGZ1 -c > $TMP.1
gunzip $PROTGZ2 -c > $TMP.2


function run
# Time: 1 min.
{
  QUERY=$1
  DB=$2
  OUT=$3
  makeblastdb  -in $DB  -dbtype prot  -blastdb_version 4  -logfile /dev/null
  blastp  -db $DB  -query $QUERY   -seg no  -comp_based_stats 0  -outfmt '6 qseqid sseqid positive nident' | sort  -k 3,3 -n -r  -k 4,4 -n -r > $OUT
}


run $TMP.1 $TMP.2 $TMP.1-2
run $TMP.2 $TMP.1 $TMP.2-1
N1=`grep -c "^>" $TMP.1`
N2=`grep -c "^>" $TMP.2`
$THIS/symbet_blastp $TMP.1-2 $TMP.2-1 $N1 $N2


rm $TMP*
