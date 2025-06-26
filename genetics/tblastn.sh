#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 7 ]; then
  echo "Run tblastn"
  echo "#1: query protein FASTA file"
  echo "#2: subject nucleotide FASTA file, may be gzip'ped"
  echo "#3: tblastn search parameters (except: -query -subject -db -out -outfmt)"
  echo "#4: tblastn output parameters (-outfmt 6 ...)"
  echo "#5: part_size"
  echo "#6: cores"
  echo "#7: output file"
  exit 1
fi
PROT=$1
DNA=$2
SEARCH=$3
OUTFMT=$4
PART=$5
CORES=$6
OUT=$7


TMP=$( mktemp )
comment $TMP


if [ -z ${DNA##*.gz} ]; then
  gunzip -c $DNA > $TMP
  DNA=$TMP
fi

mkdir $TMP.dna
# PAR
$THIS/splitFasta $DNA $TMP.dna  -group_size $PART  -mono_nuc_max 20  
mkdir $TMP.out
$THIS/../trav $TMP.dna  -step 1  -threads $CORES  "tblastn  -query $PROT  -subject %d/%f  $SEARCH  -outfmt '6 $OUTFMT' > $TMP.out/%f" 
#tblastn  -query $PROT  -subject $DNA  -comp_based_stats 0  -evalue 1e-10  -seg no  -max_target_seqs 10000  -task tblastn-fast  -threshold 100  -window_size 15  -db_gencode 11  -outfmt '6 sseqid qseqid sstart send slen qstart qend qlen nident length'
#                                                                                                                                                                                           1      2      3      4    5    6      7    8    9      10
$THIS/../trav $TMP.out "cat %d/%f" > $OUT


rm -r $TMP*
