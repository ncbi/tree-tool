#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../cpp/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Run tblastn"
  echo "#1: query protein FASTA file"
  echo "#2: subject nucleotide FASTA file, may be gzip'ped"
  exit 1
fi
PROT=$1
DNA=$2


TMP=$( mktemp )
#comment $TMP


if [ -z ${DNA##*.gz} ]; then
  gunzip -c $DNA > $TMP
  DNA=$TMP
fi
tblastn  -query $PROT  -subject $DNA  -comp_based_stats 0  -evalue 1e-10  -seg no  -max_target_seqs 10000  -task tblastn-fast  -threshold 100  -window_size 15  -db_gencode 11  -outfmt '6 sseqid qseqid sstart send slen qstart qend qlen nident length'
#                                                                                                                                                                                          1      2      3      4    5    6      7    8    9      10


rm $TMP*
