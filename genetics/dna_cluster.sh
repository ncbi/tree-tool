#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Cluster DNA sequences, if #1 = #2 then print representatives else print #1-#2 pairs"
  echo "Output: BLAST database of #2"
  echo "Temporary: #1.dir, #1.pairs"
  echo "#1; query   DNA FASTA"
  echo "#2; subject DNA FASTA"
  echo "#3: min. fraction of identity"
  echo "#4: min. fraction of coverage"
  echo "#5: output file with all pairs of sequences of the same clusters"
  echo "#6: output file with clusters: <item> <cluster representative>"
  exit 1
fi
QUERY=$1
SUBJ=$2
IDENT=$3
COV=$4
PAIR=$5
CLUST=$6


TMP=`mktemp`
echo $TMP > /dev/stderr 
#set -x


makeblastdb  -in $SUBJ   -dbtype nucl  -blastdb_version 4  -logfile /dev/null

mkdir $QUERY.dir
$THIS/splitFasta $QUERY $QUERY.dir  -group 100 

mkdir $QUERY.pairs
$THIS/../grid_wait.sh 1
$THIS/../trav $QUERY.dir "$QSUB_5 -N j%f %Q$THIS/dna_similar.sh %d/%f $SUBJ $IDENT $COV 4 > $QUERY.pairs/%f%Q > /dev/null"
$THIS/../qstat_wait.sh 3600 1

rm -r $QUERY.dir/ &

# $PAIR
if [ $QUERY == $SUBJ ]; then
  $THIS/fa2list.sh $QUERY | sed 's/^\(.*\)$/\1\t\1/1' > $PAIR
fi
$THIS/../trav $QUERY.pairs "cat %d/%f" >> $PAIR

rm -r $QUERY.pairs/ &

$THIS/../connectPairs $PAIR $CLUST  -pairs 


rm $TMP*
wait
