#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# != 3 ]; then
  echo "Find partial or extended proteins"
  echo "#1: Protein file in FASTA format"
  echo "#2: Output: Partial or extended proteins"
  echo "#3: Output: Unresolved proteins"
  exit 1
fi
IN=$1
PART_EXT=$2
UNRES=$3


TMP=$( mktemp )
#comment $TMP


cp $IN $TMP.fa

makeblastdb  -in $TMP.fa  -dbtype prot  -logfile /dev/null

# 6 min. for 4000 proteins
section "blastp"
blastp  -task blastp-fast \
  -db $TMP.fa \
  -query $TMP.fa \
  -show_gis -word_size 6 -threshold 21 -evalue 1e-20 -comp_based_stats 0 -num_threads 8 \
  -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen' \
  > $TMP.blastp

cat $TMP.blastp | awk '$4 > 60 && $4/$3 > 0.95 && $5 > $8 && $8 == 1' > $TMP.partial


section "Processing"
echo -e "#name\tdecision\tcount" > $PART_EXT
echo -e "#qseqid\tsseqid\tlength\tnident\tqstart\tqend\tqlen\tsstart\tsend\tslen" > $UNRES
while [ -s $TMP.partial ]; do
  wc -l $TMP.partial  

	cut -f 1 $TMP.partial > $TMP.1	
	cut -f 2 $TMP.partial > $TMP.2	
	cat $TMP.1 $TMP.2 > $TMP.name
	
	sort $TMP.name | uniq -c | sort -nr | head -1 > $TMP.top
	L=( $( cat $TMP.top ) )
	COUNT=${L[0]}
	NAME="${L[1]}"
	
	if [ $COUNT == 1 ]; then
	  cat $TMP.partial >> $UNRES
	  break
	fi
	
	N_EXTD=$( grep -cw "$NAME" $TMP.1 ) || true
	N_PART=$( grep -cw "$NAME" $TMP.2 ) || true
	N_EXTD_2=$(( N_EXTD * 2 ))  # PAR
	N_PART_2=$(( N_PART * 2 ))  # PAR
	DECISION=""
	if [ $N_EXTD -gt $N_PART_2 ]; then  
	  DECISION="extended"
	fi
	if [ $N_PART -gt $N_EXTD_2 ]; then
	  DECISION="partial"
	fi
	if [ "$DECISION" ]; then
	  echo -e "$NAME\t$DECISION\t$COUNT" >> $PART_EXT
	else
		grep -w "$NAME" $TMP.partial >> $UNRES
	fi

	grep -vw "$NAME" $TMP.partial > $TMP.partial1 || true
	mv $TMP.partial1 $TMP.partial
done


rm $TMP*

