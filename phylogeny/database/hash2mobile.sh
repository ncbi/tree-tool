#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Find contigs with identical protein hash codes in different GenBank bacterial assemblies"
  echo "#1: genome/ (large directory)"
  echo "#2: Assembly 1 multi-FASTA"
  echo "#3: Assembly 2 multi-FASTA"
  echo "#4: Output list of contigs in #1 common with #2"
  echo "#5: Output list of contigs in #2 common with #1"
  exit 1
fi
GENOME_DIR=$1
ASM1=$2
ASM2=$3
OUT1=$4
OUT2=$5


TMP=`mktemp`
#echo $TMP 


NAME1=`basename $ASM1`
NAME2=`basename $ASM2`

HASH1=`$THIS/../../file2hash $NAME1`
HASH2=`$THIS/../../file2hash $NAME2`

$THIS/../../setIntersect.sh $GENOME_DIR/$HASH1/$NAME1/$NAME1.hash-PRT $GENOME_DIR/$HASH2/$NAME2/$NAME2.hash-PRT 1 > $TMP.hash_common

prodigal  -o /dev/null  -i $ASM1  -d $TMP.cds1 &> /dev/null
prodigal  -o /dev/null  -i $ASM2  -d $TMP.cds2 &> /dev/null

$THIS/../../genetics/fasta2hash  $TMP.cds1  $TMP.out  -gene_finder "prodigal"  -cds  -translate  -target_hashes $TMP.hash_common  -target_seq $TMP.contig1 &> /dev/null
$THIS/../../genetics/fasta2hash  $TMP.cds2  $TMP.out  -gene_finder "prodigal"  -cds  -translate  -target_hashes $TMP.hash_common  -target_seq $TMP.contig2 &> /dev/null

cut -f 1 -d ' ' $TMP.contig1 | sed 's/_[0-9]\+$//1' | sort -u > $OUT1
cut -f 1 -d ' ' $TMP.contig2 | sed 's/_[0-9]\+$//1' | sort -u > $OUT2


rm -r $TMP* 
