#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find proteins with introns by tblastn"
  echo "#1: DNA multi-FASTA file"
  echo "#2: reference proteins"
  echo "#3: output protein file"
  echo "Time: 110 min./3.8G"
  exit 1
fi
DNA=$1
REF=$2
PROT=$3


TMP=$( mktemp )
#comment $TMP  
#set -x  

MATRIX="BLOSUM62"

makeblastdb  -in $DNA  -out $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
tblastn  -query $REF  -db $TMP  -db_gencode 1  -seg no  -comp_based_stats 0  -word_size 3  -evalue 1000  -matrix $MATRIX  -outfmt '6 qseqid sseqid qstart qend sstart send qseq sseq' -num_threads 4 > $TMP.tblastn
  # >12 hours:
    # -word_size 2
    # -window_size 0
    # -threshold 10
  # From Terence M. (Eukaryotic annotation):
  #   blastp_wnode  -best_hit_overhang 0.1  -best_hit_score_edge 0.1  -comp_based_stats 0  -evalue 0.0001  -seg no  -threshold 21  -word_size 6  \
  #                 -asn-cache ${GP_cache_dir}  -backlog 1  -delay 0  -max-jobs 1  -service ${GP_qservice}  -workers 4
  #   diamond  --comp-based-stats 0  --evalue 0.0001  --very-sensitive  --masking 0  --unal 0  seq-align-set  \
  #            -asn-cache ${GP_cache_dir}  -blastp-args  --sam-query-len  -diamond-executable ${GP_HOME}/third-party/diamond/diamond  -lds2 ${input.lds2.target}  -ofmt  -output-dir ${output}  -output-manifest ${output}/align.mft -output-prefix hits -query-fmt seq-ids -query-manifest ${input.query_ids} -subject-fmt seq-ids -subject-manifest ${input.subject_ids} -work-area ${tmp}

$THIS/tblastn2marker_euk  $TMP.tblastn  -matrix $MATRIX  -qc  -log $TMP.log > $PROT
# Quality: 
  # grep '^>' $PROT | sed 's/^.*score=//1' | count
  # fasta2len $PROT | cut -f 2 | count


rm $TMP*  
