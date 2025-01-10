#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Run tblastn"
  echo "#1: query protein FASTA file"
  echo "#2: subject nucleotide FASTA file, may be gzip'ped"
  echo "#3: tblastn parameters (except: -query -subject -db -out)"
  echo "#4: cores"
  echo "#5: output file"
  exit 1
fi
PROT=$1
DNA=$2
PARM=$3
CORES=$4
OUT=$5


TMP=$( mktemp )
comment $TMP


if [ -z ${DNA##*.gz} ]; then
  gunzip -c $DNA > $TMP
  DNA=$TMP
fi
mkdir $TMP.dna
$THIS/splitFasta $DNA $TMP.dna -group_size 50000000
mkdir $TMP.out
$THIS/../trav $TMP.dna  -step 1  -threads $CORES  "tblastn  -query $PROT  -subject %d/%f  $PARM  > $TMP.out/%f" 
#tblastn  -query $PROT  -subject $DNA  -comp_based_stats 0  -evalue 1e-10  -seg no  -max_target_seqs 10000  -task tblastn-fast  -threshold 100  -window_size 15  -db_gencode 11  -outfmt '6 sseqid qseqid sstart send slen qstart qend qlen nident length'
#                                                                                                                                                                                           1      2      3      4    5    6      7    8    9      10
$THIS/../trav $TMP.out "cat %d/%f" > $OUT


rm -r $TMP*
