#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# != 3 ]; then
  echo "MLST-hash typing"
  echo "Print result"
  echo "#1: DNA identifier"
  echo "#2: DNA sequence BLAST database"
  echo "#3: MLST FASTA file with different loci"
  exit 1
fi
ID=$1
DNA=$2
MLST=$3


if [ ! $DNA ]; then
  error "Empty DNA file name"
fi

if [ ! -e $DNA ]; then
  error "File $DNA does not exist"
fi

if [ ! -s $DNA ]; then
  error "File $DNA is empty"
fi


TMP=$( mktemp )
#comment $TMP   


N=$( grep -c '^>' $MLST )
# PAR
blastn  -query $MLST  -db $DNA  -dust no  -evalue 1e-50  -outfmt '6 qseqid sseq nident length' | awk '$3/$4 > 0.8' | awk '{print $1, $2;}' > $TMP.blastn
$THIS/mlst2hash $N $TMP.blastn $TMP.type
echo "$ID" $( cat $TMP.type )


rm -f $TMP*  
