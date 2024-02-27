#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../cpp/bash_common.sh
if [ $# != 3 ]; then
  echo "Find partial or extended proteins"
  echo "#1: Protein file in FASTA format"
  echo "#2: Output: Partial or extended proteins"
  echo "#3: Output: Unresolved proteins"
  exit 1
fi


TMP=`mktemp`


# Make all start codons to be 'M'
cat $1 | sed '/^>/\!s/[LIV]/M/g' > $TMP.fa

makeblastdb  -in $TMP.fa  -dbtype prot  -logfile /dev/null

# 6 min. for 4000 proteins
section "Blastp"
blastp  -task blastp-fast \
  -db $TMP.fa \
  -query $TMP.fa \
  -show_gis -word_size 6 -threshold 21 -evalue 1e-20 -comp_based_stats 0 -num_threads 8 \
  -outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen' \
  > $TMP.blastp

cat $TMP.blastp | awk '$4 > 60 && $4/$3 > 0.95 && $5 > $8 && $8 == 1' > $TMP.partial


section "Processing"
echo -e "#name\tdecision\tcount" > $2
echo -e "#qseqid\tsseqid\tlength\tnident\tqstart\tqend\tqlen\tsstart\tsend\tslen" > $3
while [ -s $TMP.partial ]; do
  wc -l $TMP.partial  

	cut -f 1 $TMP.partial > $TMP.1	
	cut -f 2 $TMP.partial > $TMP.2	
	cat $TMP.1 $TMP.2 > $TMP.name
	
	sort $TMP.name | uniq -c | sort -nr | head -1 > $TMP.top
	L=(`cat $TMP.top`)
	COUNT=${L[0]}
	NAME="${L[1]}"
	
	if [ $COUNT == 1 ]; then
	  cat $TMP.partial >> $3
	  break
	fi
	
	N_EXTD=`grep -cw "$NAME" $TMP.1`
	N_PART=`grep -cw "$NAME" $TMP.2`
	N_EXTD_2 = $(( $N_EXTD * 2 ))  # PAR
	N_PART_2 = $(( $N_PART * 2 ))   # PAR
	DECISION=""
	if [ $N_EXTD -gt $N_PART_2) ]; then  
	  DECISION="extended"
	fi
	if [ $N_PART -gt $N_EXTD_2 ]; then
	  DECISION="partial"
	fi
	if [ "$DECISION" ]; then
	  echo "$NAME $DECISION $COUNT" >> $2
	else
		grep -w "$NAME" $TMP.partial >> $3
	endif

	grep -vw "$NAME" $TMP.partial > $TMP.partial1
	mv $TMP.partial1 $TMP.partial
done


rm -f $TMP*

