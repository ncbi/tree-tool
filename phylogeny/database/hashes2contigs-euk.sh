#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Find contigs with identical protein hash codes in two GenBank eukaryotic assemblies"
  echo "#1: genome/ (large directory)"
  echo "#2: assembly 1 FASTA (named as fullasm_id)"
  echo "#3: assembly 2 FASTA (named as fullasm_id)"
  echo "#4: output #2-#3.ctg"
  echo "#5: output #3-#2.ctg"
  exit 1
fi
GENOME_DIR=$1
ASM=($2 $3)
OUT=($4 $5)


TMP=$( mktemp )
comment $TMP 
#set -x


NAME=()
HASH=()
for i in {0..1}; do
  NAME+=( $( basename ${ASM[$i]} ) )
  HASH+=( $( $THIS/../../file2hash ${NAME[$i]} ) )
done

$THIS/../../setIntersect.sh $GENOME_DIR/${HASH[0]}/${NAME[0]}/${NAME[0]}.hash-PRT $GENOME_DIR/${HASH[1]}/${NAME[1]}/${NAME[1]}.hash-PRT 1 > $TMP.hash_common
wc -l $TMP.hash_common


function run 
{
  local i=$1
  section "Making ${OUT[$i]}"
  $THIS/../../check_file.sh $GENOME_DIR/${HASH[$i]}/${NAME[$i]}/genemark.gtf.gz 1
  zcat $GENOME_DIR/${HASH[$i]}/${NAME[$i]}/genemark.gtf.gz > $TMP.genemark.gtf.$i
  $THIS/../../genetics/GeneMark2CDS ${ASM[$i]} $TMP.genemark.gtf.$i  -gtf  -gencode 1  -cds $TMP.cds.$i  -prot_len_min 150  -complete  -qc &> /dev/null
  $THIS/../../genetics/fasta2hash  $TMP.cds.$i  $TMP.out  -gene_finder "GeneMark"  -cds  -translate  -target_hashes $TMP.hash_common  -target_seq $TMP.contig.$i  &> /dev/null
  echo -e "#CDSs\tContig" > $TMP.contig.tsv
  cut -f 1 -d ' ' $TMP.contig.$i | sed 's/\.[0-9]\+$//1' | sort | uniq -c | sed 's/^ *//1' | tr ' ' '\t' >> $TMP.contig.tsv
  echo -e "#Contig\tlen" > $TMP.out
  $THIS/../../genetics/fasta2len ${ASM[$i]} 1>> $TMP.out 2> /dev/null
  $THIS/../../tsv/tsv_expand.sh $TMP.out $TMP.contig.tsv '' &> /dev/null
  mv $TMP.out ${OUT[$i]}
}

run 0
run 1


rm -r $TMP* 
