#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 8 ]; then
  echo "Create annotation files for #1 in the current directory"
  echo "Invokes GeneMark"
  echo "#1: eukaryotic genome name"
  echo "#2: eukaryotic DNA FASTA"
  echo "#3: kingdom: Fungi|Viridiplantae|''"
  echo "#4: universal HMM library (absolute path) or ''"
  echo "#5: Pfam HMM library (absolute path) or ''"
  echo "#6: use Pfam HMM cutoff (0/1)"
  echo "#7: number of cores"
  echo "#8: log file (absolute path)"
  echo "Time: #1 size: 10M -1 hour; 4G - 6 days"
  exit 1
fi
ASM=$1
FASTA=$2
KINGDOM="$3"
UNIV=$4
PFAM=$5
PFAM_CUTOFF=$6
CORES=$7
LOG=$8


#set -x

if [ $KINGDOM -a $KINGDOM != "Fungi" -a $KINGDOM != "Viridiplantae" ]; then
  error "Unknown kingdom '$KINGDOM'"
fi


if [ ! -e $FASTA ]; then
  error "No $FASTA" >> $LOG
fi
if [ ! -s $FASTA ]; then
  error "Empty genome" >> $LOG
fi

DIR=`dirname $FASTA`

cd $DIR


ls -laF $FASTA
$THIS/dna2stat $FASTA  -log $LOG > $ASM.stat


if [ ! -e $ASM.prot ]; then
  echo "GeneMark-ES 4.69" > annot_software

  # genemark.gtf
  FUNGUS_PAR=""
  if [ $KINGDOM == "Fungi" ]; then
    FUNGUS_PAR="--fungus"
  fi
  $THIS/../rm_all.sh data
  $THIS/../rm_all.sh info
  $THIS/../rm_all.sh output
  $THIS/../rm_all.sh run
  rm -f gmes.log
  rm -f run.cfg
  # Approximate
  /usr/local/gmes/4.69/bin/gmes_petap.pl  --ES  $FUNGUS_PAR  --sequence $FASTA  --soft_mask 0  --cores $CORES  &>> $LOG
  echo -e "\nAnnotation finihsed!\n" >> $LOG
  $THIS/../rm_all.sh data
  $THIS/../rm_all.sh info
  $THIS/../rm_all.sh output
  $THIS/../rm_all.sh run
  rm gmes.log
  rm run.cfg
  rm gmhmm.mod

  # PAR
  PROT_LEN=150

  section "$ASM.prot"
  $THIS/GeneMark2CDS $FASTA genemark.gtf  -gtf  -gencode 1  -prot $ASM.prot  -prot_len_min 20  -ambig 10  -log $LOG  -qc
    # was: -prot_len_min 60

  section "$ASM.cds"
  $THIS/GeneMark2CDS $FASTA genemark.gtf  -gtf  -gencode 1  -cds $ASM.cds  -prot_len_min $PROT_LEN  -complete  -log $LOG  -qc
  gzip genemark.gtf

  section "$ASM.hash-CDS"
  $THIS/fasta2hash $ASM.cds $ASM.hash-CDS  -cds  -gene_finder "GeneMark"  -prot_len_min $PROT_LEN  -log $LOG  -qc

  section "$ASM.hash-PRT"
  $THIS/fasta2hash $ASM.cds $ASM.hash-PRT  -cds  -gene_finder "GeneMark"  -prot_len_min $PROT_LEN  -log $LOG  -translate  -qc  
  rm $ASM.cds
fi


if [ $PFAM ]; then
  section "Pfam"
  $THIS/prots2hmm_hash.sh $ASM.prot $PFAM $PFAM_CUTOFF $ASM.HMM $ASM.hash-HMM $CORES $LOG >> $LOG
  gzip $ASM.HMM
fi

if [ $UNIV -a $KINGDOM != "Viridiplantae" ]; then
  section "prots2hmm_univ.sh"
  $THIS/prots2hmm_univ.sh $ASM $UNIV 0 $CORES $LOG >> $LOG
fi


gzip $ASM.prot


rm -f $LOG


