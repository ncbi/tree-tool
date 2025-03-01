#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# != 9 ]; then
  echo "Print Genome.id's approximately closest to #1, where Genome.in_tree = 1"
  echo "#1: SQL Server name"
  echo "#2: database with tables Genome, GenomeHash, FreqHash"
  echo "#3: directory for bulk insert"
  echo "#4: path in Universal Naming Convention to the bulk insert directory"
  echo "#5: Genome.id"
  echo "#6: Taxroot.id"
  echo "#7: directory with #1/#1.hash-{[CDS,]PRT,HMM}"
  echo "#8: Hashtype.id for lowest-level dissimilarity"
  echo "#9: subset of Genome.id for search"
  exit 1
fi
SERVER=$1
DATABASE=$2
BULK=$3
BULK_REMOTE=$4
GENOME=$5
TAXROOT=$6
DIR=$7
TYPE=$8
SUBSET=$9


TMP=$( mktemp )
#comment $TMP
#set -x


NAME=$( basename $TMP )
BULK_FILE=$BULK/$NAME
FILE=$DIR/$GENOME

cp $SUBSET $BULK_FILE.subset
unix2dos -o $BULK_FILE.subset &> /dev/null

cp $FILE.hash-$TYPE $BULK_FILE
unix2dos -o $BULK_FILE &> /dev/null
$THIS/Genome_hash_requestClosest_.sh $SERVER $DATABASE $BULK_REMOTE $GENOME $TAXROOT $TYPE $NAME $NAME.subset > $TMP

if [ "$TYPE" != "PRT" ]; then
  N=$( < $TMP  wc -l )
  if [ $N -lt 90 ]; then  # PAR
    cp $FILE.hash-PRT $BULK_FILE
    unix2dos -o $BULK_FILE &> /dev/null
    $THIS/Genome_hash_requestClosest_.sh $SERVER $DATABASE $BULK_REMOTE $GENOME $TAXROOT 'PRT' $NAME $NAME.subset >> $TMP
		sort -u $TMP > $TMP.1
		mv $TMP.1 $TMP
  fi
fi

N=$( < $TMP  wc -l)
if [ $N -lt 90 ]; then  # PAR
  cp $FILE.hash-HMM $BULK_FILE
  unix2dos -o $BULK_FILE &> /dev/null
  $THIS/Genome_hash_requestClosest_.sh $SERVER $DATABASE $BULK_REMOTE $GENOME $TAXROOT 'HMM' $NAME $NAME.subset >> $TMP
fi

sort -u $TMP 


rm $BULK_FILE*
rm $TMP*
