#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Find proteins with introns by tblastn"
  echo "#1: DNA multi-FASTA file"
  echo "#2: reference proteins"
  echo "#3: output protein file"
  exit 1
fi
DNA=$1
REF=$2
PROT=$3


TMP=$( mktemp )
comment $TMP  
#set -x  

makeblastdb  -in $DNA  -out $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
tblastn  -query $REF  -db $TMP  -db_gencode 1   -outfmt '6 qseqid sseqid qstart qend sstart send qseq sseq' -num_threads 4 > $TMP.tblastn
  # -seg no  -comp_based_stats 0 
  # (!) -evalue 1000  -word_size 2 

# From Terence M.:
# blastp_wnode -asn-cache ${GP_cache_dir} -backlog 1 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -comp_based_stats 0 -delay 0 -evalue 0.0001 -max-jobs 1 -seg no -service ${GP_qservice} -threshold 21 -word_size 6 -workers 4
# iamond -asn-cache ${GP_cache_dir} -blastp-args --sam-query-len --comp-based-stats 0 --evalue 0.0001 --very-sensitive --masking 0 --unal 0 -diamond-executable ${GP_HOME}/third-party/diamond/diamond -lds2 ${input.lds2.target} -ofmt seq-align-set -output-dir ${output} -output-manifest ${output}/align.mft -output-prefix hits -query-fmt seq-ids -query-manifest ${input.query_ids} -subject-fmt seq-ids -subject-manifest ${input.subject_ids} -work-area ${tmp}

# Nick: -word_size 3  -evalue 100

$THIS/tblastn2marker_euk  $TMP.tblastn  -qc  -log $TMP.log > $PROT


rm $TMP*  
