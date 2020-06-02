#!/bin/bash 
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Create files for #1 in the current directory"
  echo "Invokes: prodigal"
  echo "#1: prokaryotic genome name"
  echo "#2: DNA FASTA of a prokaryotic genome"
  echo "#3: HMM_CUTOFF (0/1/null; null - do not use HMMs)"
  echo "#4: MLST scheme taxid | 0"
  echo "#5: log file"
  echo "#6: Keep proteins (0/1)"
  echo "Time: 6 min."
  exit 1
fi
NAME=$1
FASTA=$2
HMM_CUTOFF=$3
MLST=$4
LOG=$5
KEEP_PROT=$6


GeneMark=0  # PAR


$THIS/dna2stat $FASTA  -log $LOG > $NAME.stat


# GENE_FINDER, $NAME.cds
if [ $GeneMark == 1 ]; then
	GENE_FINDER=GeneMark
  # GeneMarkS-2
	/home/brovervv/gms2_v102_ncbi/gms2.pl  --seq $FASTA  --genome-type bacteria  --mgm-type bac  --fnn $NAME.cds  &>> $LOG
	rm log
	rm *.lst
	rm *.mod
else
	GENE_FINDER="prodigal"
	prodigal  -o /dev/null  -i $FASTA  -d $NAME.cds &>> $LOG
fi

$THIS/fasta2hash $NAME.cds $NAME.hash-CDS  -log $LOG  -cds  -gene_finder $GENE_FINDER              -min_prot_len 150  
$THIS/fasta2hash $NAME.cds $NAME.hash-PRT  -log $LOG  -cds  -gene_finder $GENE_FINDER  -translate  -min_prot_len 150  

if [ $HMM_CUTOFF != "null" ]; then
  $THIS/fasta2hash $NAME.cds $NAME.hash-PRT1 -log $LOG  -cds  -gene_finder $GENE_FINDER  -translate  -min_prot_len 0  -out_prot $NAME.prot
  rm $NAME.hash-PRT1
  $THIS/prots2hmm_univ.sh $NAME /home/brovervv/panfs/GenBank/bacteria/hmm-univ.LIB 1 $LOG
  echo "Pfam ..."
  $THIS/prots2hmm_hash.sh $NAME.prot /home/brovervv/panfs/Pfam/Pfam-A.hmm $HMM_CUTOFF $NAME.HMM $NAME.hash-HMM $LOG
  gzip $NAME.HMM
  if [ $KEEP_PROT == 1 ]; then
    gzip $NAME.prot
  else
    rm $NAME.prot
  fi
fi
rm $NAME.cds


if [ $MLST == 0 ]; then
  cp /dev/null $NAME.mlst
else
  $THIS/genome2mlst.sh $NAME $FASTA ~brovervv/panfs/MLST/$MLST > $NAME.mlst
fi


