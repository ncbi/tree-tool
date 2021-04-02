#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# != 7 ]; then
  echo "Print Genome.id's approximately closest to #1, where Genome.in_tree = 1"
  echo "#1: SQL Server name"
  echo "#2: Database with tables Genome, GenomeHash, FreqHash"
  echo "#3: Directory for bulk insert"
  echo "#4: path in Universal Naming Convention to the bulk insert directory"
  echo "#5: Genome.id"
  echo "#6: directory with #1/#1.hash-{[CDS,]PRT,HMM}"
  echo "#7: Hashtype.id for lowest-level dissimilarity"
  exit 1
fi
SERVER=$1
DATABASE=$2
BULK=$3
BULK_REMOTE=$4
GENOME=$5
DIR=$6
TYPE=$7


TMP=`mktemp`  
#echo $TMP
#set -x


NAME=`basename $TMP`
BULK_FILE=$BULK/$NAME

FILE=$DIR/$GENOME/$GENOME

cp $FILE.hash-$TYPE $BULK_FILE
unix2dos -o $BULK_FILE &> /dev/null
$THIS/Genome_hash_requestClosest_.sh $SERVER $DATABASE $BULK_REMOTE $GENOME $TYPE $NAME > $TMP

if [ "$TYPE" != "PRT" ]; then
  N=`cat $TMP | wc -l`
  if [ $N -lt 90 ]; then  # PAR
    cp $FILE.hash-PRT $BULK_FILE
    unix2dos -o $BULK_FILE &> /dev/null
    $THIS/Genome_hash_requestClosest_.sh $SERVER $DATABASE $BULK_REMOTE $GENOME 'PRT' $NAME >> $TMP
		sort -u $TMP > $TMP.1
		mv $TMP.1 $TMP
  fi
fi

N=`cat $TMP | wc -l`
if [ $N -lt 90 ]; then  # PAR
  cp $FILE.hash-HMM $BULK_FILE
  unix2dos -o $BULK_FILE &> /dev/null
  $THIS/Genome_hash_requestClosest_.sh $SERVER $DATABASE $BULK_REMOTE $GENOME 'HMM' $NAME >> $TMP
fi

sort -u $TMP 


rm $TMP*
rm $BULK_FILE
