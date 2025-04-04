#!/bin/bash --noprofile
# PD-2992
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
HEADER="#Target\tBlastRule\tBlastRule_protein"
if [ $# -ne 3 ]; then
  # PD-2992
  echo "Find Blast Rules hitting the complete input proteins"
  echo "Print: $HEADER, sorted by Target"
  echo "#1: input protein(s) file"
  echo "#2: Blast Rule data directory"
  echo "#3: output report"
  exit 1
fi
IN=$1
DB=$2
OUT=$3


TMP=$( mktemp )
#comment $TMP  


blastp  -query $IN  -db $DB/prot  -show_gis  -comp_based_stats 0  -evalue 1e-10  -num_threads 10 \
   -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq sseq' \
   -out $TMP.blastp

sort -k 1,2 $TMP.blastp > $TMP.blastp-sorted
# PAR
$THIS/blastp_merge $TMP.blastp-sorted 0.5 20 | sort -k 2 > $TMP.merge
# qseqid sseqid ident_frac p_target_coverage p_ref_coverage"

join  -1 2  -2 1  $TMP.merge $DB/prot_br | sort -k 6 > $TMP.blastp_br

join  -1 6  -2 1  $TMP.blastp_br $DB/br > $TMP.blastp_br_cutoff
# NBR000415 1215061635 ACB17037.1 20.97 44.46      64.33 BlastRuleException 1             94          90                    90                    96                 90                   25
#      1        2          3       4     5          6            7          8              9          10                    11                    12                 13                   14
#     BR       ref      target  ident target_cov ref_cov      type        precedence complete_ident complete_wp_coverage  complete_br_coverage  partial_ident  partial_wp_coverage   partial_br_coverage

cat $TMP.blastp_br_cutoff | awk '$4 >= $9 && $5 >= $10 && $6 >= $11' | cut -f 1,2,3,8 -d ' ' | sort -k 4 -n | awk '{print $3, $1, $2, $4};' | sort -k 1 > $TMP.pass


echo -e "$HEADER" > $OUT
TARGET_PREV=""
PRECEDENCE_PREV=0
while read -r LINE
do
  ARR=( $LINE )
  
  TARGET=${ARR[0]}
  BR=${ARR[1]}
  REF=${ARR[2]}
  PRECEDENCE=${ARR[3]}
  
  if [ "$TARGET" \< "$TARGET_PREV" ]; then
    error "Target is not sorted"
  fi
  
  if [ "$TARGET" == "$TARGET_PREV" ]; then
    if [ $PRECEDENCE == $PRECEDENCE_PREV ]; then
      echo -e "$TARGET\t$BR\t$REF" >> $OUT
    fi
  else
    echo -e "$TARGET\t$BR\t$REF" >> $OUT
  fi
  
  TARGET_PREV=$TARGET
  PRECEDENCE_PREV=$PRECEDENCE
done < $TMP.pass


rm $TMP* 

