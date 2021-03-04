#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print <start> <end> <strand> of protein including stop codon in #2 if matched"
  echo "#1: protein sequence file (1 line, not FASTA, without stop codon)"
  echo "#2: DNA FASTA file"
  exit 1
fi
PROT=$1
DNA=$2


TMP=`mktemp`
#echo $TMP


echo '>seq' >  $TMP.prot
cat $PROT   >> $TMP.prot
echo '*'    >> $TMP.prot

cp $DNA $TMP.dna
makeblastdb  -in $TMP.dna  -dbtype nucl  -logfile /dev/null  -blastdb_version 4

tblastn  -query $TMP.prot   -db $TMP.dna  -word_size 3  -seg no  -db_gencode 11  -outfmt '6 qstart qend qlen sstart send length nident sseq' 2> /dev/null | sort -k 4,4 -n | sort -k 7,7 -n -r | head -1 > $TMP.tblastn 
RES=(`cat $TMP.tblastn`)
QSTART=${RES[0]}
QEND=${RES[1]}
QLEN=${RES[2]}
SSTART=${RES[3]}
SEND=${RES[4]}
LENGTH=${RES[5]}
NIDENT=${RES[6]}
SSEQ="${RES[7]}"

STRAND=1
if [ $SSTART -gt $SEND ]; then
  STRAND=0
  SSTART_=$SSTART
  SSTART=$SEND
  SEND=$SSTART_
fi

if [ $QSTART == 1 ]; then
  FIRST_CHAR=`echo $SSEQ | cut -c1`
  if [ $FIRST_CHAR == 'L' -o $FIRST_CHAR == 'I' -o $FIRST_CHAR == 'V' ]; then
    NIDENT=$(( $NIDENT + 1 ))
  fi
fi

if [ $LENGTH == $NIDENT -a $LENGTH == $QLEN -a $QSTART == 1 -a $QEND == $QLEN ]; then
  echo "$SSTART $SEND $STRAND"  
fi
  

rm $TMP*
