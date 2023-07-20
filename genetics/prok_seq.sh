#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Create files for #1 in the current directory"
  echo "Invokes: prodigal"
  echo "#1: prokaryotic genome name"
  echo "#2: DNA FASTA of a prokaryotic genome"
  echo "#3: universal HMM library or ''"
  echo "#4: Pfam HMM library  or ''"
  echo "#5: use Pfam HMM cutoff (0/1)"
  echo "#6: log file"
  echo "Time: 6 min."
  exit 1
fi
NAME=$1
FASTA=$2
UNIV=$3
PFAM=$4
PFAM_CUTOFF=$5
LOG=$6


$THIS/dna2stat $FASTA  -log $LOG > $NAME.stat


# GENE_FINDER, $NAME.cds
if false; then
	GENE_FINDER="GeneMark"
  # GeneMarkS-2
	$HOME/gms2_v102_ncbi/gms2.pl  --seq $FASTA  --genome-type bacteria  --mgm-type bac  --fnn $NAME.cds  &>> $LOG
	rm log
	rm *.lst
	rm *.mod
else
	GENE_FINDER="prodigal"
	prodigal  -o /dev/null  -i $FASTA  -d $NAME.cds  &>> $LOG
fi

MIN_PROT_LEN=150
$THIS/fasta2hash $NAME.cds $NAME.hash-CDS  -log $LOG  -cds  -gene_finder $GENE_FINDER              -prot_len_min $MIN_PROT_LEN  
$THIS/fasta2hash $NAME.cds $NAME.hash-PRT  -log $LOG  -cds  -gene_finder $GENE_FINDER  -translate  -prot_len_min $MIN_PROT_LEN  

$THIS/fasta2hash $NAME.cds $NAME.hash-PRT1 -log $LOG  -cds  -gene_finder $GENE_FINDER  -translate  -prot_len_min 0  -ambig 10  -out_prot $NAME.prot
rm $NAME.hash-PRT1

rm $NAME.cds

# PAR
CORES=4

if [ $UNIV ]; then
  section "univ"
  $THIS/prots2hmm_univ.sh $NAME $UNIV 1 $CORES $LOG
fi

if [ $PFAM ]; then
  section "Pfam"
  $THIS/prots2hmm_hash.sh $NAME.prot $PFAM $PFAM_CUTOFF $NAME.HMM $NAME.hash-HMM $CORES $LOG
  gzip $NAME.HMM
fi

# Find mobilome (transposases) ??

gzip $NAME.prot



