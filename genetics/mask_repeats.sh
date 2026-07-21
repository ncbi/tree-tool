#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find and mask repeats in #1"
  echo "#1: DNA multi-FASTA"
  echo "#2: threads"
  echo "#3: output file with the coordinates of masked segments"
  echo "#4: output masked #1"
  exit 1
fi
ASM=$1
THREADS=$2
MSK=$3
OUT=$4


# PAR
FRAG=1000
HANG=$(( FRAG / 10 ))


TMP=$( mktemp )
comment $TMP 


section "Fragmentation by $FRAG bp"
$THIS/../dir_hash_init.sh $TMP.frag 1000
mkdir $TMP.seq
$THIS/splitFasta $ASM $TMP.seq  -qc
ls $TMP.seq > $TMP.seqs
wc -l $TMP.seqs
N=0
while read -r F
do
  N=$(( N + 1 ))
  printf "\r%s" "$N"
  $THIS/fragmentize $TMP.seq/$F $FRAG  -qc -noprogress > $TMP.fragmented
  $THIS/splitFasta $TMP.fragmented $TMP.frag  -len_min $FRAG  -large  -qc -noprogress
done < $TMP.seqs
echo ""

section "blastn"
mkdir $TMP.out
makeblastdb  -in $ASM   -dbtype nucl  -blastdb_version 4  -logfile /dev/null  -out $TMP.db
$THIS/../trav $TMP.frag "ls %d/%f" > $TMP.frags
wc -l $TMP.frags
$THIS/../trav $TMP.frags -step 1 "$THIS/find_repeats_.sh $TMP.frag/%h/%f $TMP.db $HANG $TMP.out"  -threads $THREADS  -errors $TMP.err
if [ -s $TMP.err ]; then
  error "$TMP.err"
fi
trav $TMP.out "cat %d/%f" > $TMP.msk

section "dna_mask"
MASK_MIN=$(( FRAG - 2 * HANG ))
$THIS/dna_mask $ASM $TMP.msk  -merged_mask $MSK  -mask_min $MASK_MIN  -qc > $OUT


rm -r $TMP*
