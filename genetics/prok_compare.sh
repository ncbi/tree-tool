#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Compare two close genomes by counting common long CDSs or proteins (>= 150 aa)"
  echo "Invokes: prodigal"
  echo "#1: Genome 1 (DNA FASTA)"
  echo "#2: Genome 2 (DNA FASTA)"
  exit 1
fi


TMP=`mktemp`
#echo $TMP  


prodigal  -o /dev/null  -i $1  -d $TMP.cds1 &> /dev/null 
prodigal  -o /dev/null  -i $2  -d $TMP.cds2 &> /dev/null

$THIS/fasta2hash $TMP.cds1 $TMP.hash-CDS1  -cds  -gene_finder prodigal              -prot_len_min 150 &> /dev/null
$THIS/fasta2hash $TMP.cds1 $TMP.hash-PRT1  -cds  -gene_finder prodigal  -translate  -prot_len_min 150 &> /dev/null 

$THIS/fasta2hash $TMP.cds2 $TMP.hash-CDS2  -cds  -gene_finder prodigal              -prot_len_min 150 &> /dev/null
$THIS/fasta2hash $TMP.cds2 $TMP.hash-PRT2  -cds  -gene_finder prodigal  -translate  -prot_len_min 150 &> /dev/null

CDS1=`cat $TMP.hash-CDS1 | wc -l`
CDS2=`cat $TMP.hash-CDS2 | wc -l`
CDS=`$THIS/../setIntersect.sh $TMP.hash-CDS1 $TMP.hash-CDS2 1 | wc -l`
echo "# Long CDSs in $1: $CDS1"
echo "# Long CDSs in $2: $CDS2"
echo "# Common CDSs: $CDS"

PRT1=`cat $TMP.hash-PRT1 | wc -l`
PRT2=`cat $TMP.hash-PRT2 | wc -l`
PRT=`$THIS/../setIntersect.sh $TMP.hash-PRT1 $TMP.hash-PRT2 1 | wc -l`
echo "# Long proteins in $1: $PRT1"
echo "# Long proteins in $2: $PRT2"
echo "# Common proteins: $PRT"


rm -f $TMP*


