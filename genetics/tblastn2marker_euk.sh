#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Find proteins with introns by tblastn"
  echo "#1: DNA multi-FASTA file"
  echo "#2: reference proteins, where if #3 = 1 then each id=<gene>-<variant>, where <gene> has no dashes"
  echo "#3: delimiter between gene and variant | ''"
  echo "#4: cores"
  echo "#5: output protein file"
  echo "#6: log file"
  echo "Time: 1 day/3G"
  exit 1
fi
DNA=$1
REF=$2
DELIM=$3
CORES=$4
PROT=$5
LOG=$6


TMP=$( mktemp )
comment $TMP &> $LOG
date >> $LOG
ls -laF $DNA >> $LOG
#set -x  

MATRIX="BLOSUM62"

SEARCH="-comp_based_stats 0  -seg no  -max_target_seqs 10000  -dbsize 10000  -evalue 1  -word_size 3  -matrix $MATRIX  -db_gencode 1"
  # Cf. Hsp::blastp_slow
                       
  # Too much time:
    # -word_size 2
    # -window_size 0
    # -threshold 10
  # From Terence M. (Eukaryotic annotation):
  #   blastp_wnode  -best_hit_overhang 0.1  -best_hit_score_edge 0.1  -comp_based_stats 0  -evalue 0.0001  -seg no  -threshold 21  -word_size 6  \
  #                 -asn-cache ${GP_cache_dir}  -backlog 1  -delay 0  -max-jobs 1  -service ${GP_qservice}  -workers 4
  #   diamond  --comp-based-stats 0  --evalue 0.0001  --very-sensitive  --masking 0  --unal 0  seq-align-set  \
  #            -asn-cache ${GP_cache_dir}  -blastp-args  --sam-query-len  -diamond-executable ${GP_HOME}/third-party/diamond/diamond  -lds2 ${input.lds2.target}  -ofmt  -output-dir ${output}  -output-manifest ${output}/align.mft -output-prefix hits -query-fmt seq-ids -query-manifest ${input.query_ids} -subject-fmt seq-ids -subject-manifest ${input.subject_ids} -work-area ${tmp}
#makeblastdb  -in $DNA  -out $TMP  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
#tblastn  -query $REF  -db $TMP  -db_gencode 1  -seg no  -comp_based_stats 0  -word_size 3  -evalue 1000  -matrix $MATRIX  -outfmt '6 qseqid sseqid qstart qend sstart send qseq sseq' -num_threads $CORES  > $TMP.tblastn
$THIS/tblastn.sh $REF $DNA "$SEARCH"  "qseqid sseqid qstart qend qlen sstart send slen qseq sseq"  10000000 $CORES $TMP.tblastn
sort -k 1,2 $TMP.tblastn > $TMP.tblastn-sort

DELIM_PAR=""
if [ $DELIM ]; then
  DELIM_PAR="-delimiter $DELIM"
fi
$THIS/tblastn2marker_euk  $TMP.tblastn-sort  $DELIM_PAR  -matrix $MATRIX  -threshold 3.5  -qc  -log $LOG > $PROT
# Quality: 
  # grep '^>' $PROT | sed 's/^.*score=//1' | count
  # fasta2len $PROT -noprogress | cut -f 2 | count


rm $TMP*  
