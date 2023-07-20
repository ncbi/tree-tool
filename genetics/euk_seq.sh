#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 7 ]; then
  echo "Create annotation files for #1 in the directory of #1"
  echo "#1: eukaryotic DNA FASTA"
  echo "#2: fungus (0/1)"
  echo "#3: universal HMM library (absolute path) or ''"
  echo "#4: Pfam HMM library (absolute path) or ''"
  echo "#5: use Pfam HMM cutoff (0/1)"
  echo "#6: number of cores"
  echo "#7: log file (absolute path)"
  echo "Time: #1 size: 10M -1 hour; 4G - 6 days"
  exit 1
fi
FASTA=$1
FUNGUS=$2
UNIV=$3
PFAM=$4
PFAM_CUTOFF=$5
CORES=$6
LOG=$7


#set -x


if [ ! -e $FASTA ]; then
  error "No $FASTA" >> $LOG
fi
if [ ! -s $FASTA ]; then
  error "Empty genome" >> $LOG
fi

DIR=`dirname $FASTA`
ASM=`basename $FASTA`

cd $DIR


ls -laF $ASM
$THIS/dna2stat $ASM  -log $LOG > $ASM.stat


if [ ! -e $ASM.prot ]; then
  echo "GeneMark-ES 4.69" > annot_software

  # genemark.gtf
  FUNGUS_PAR=""
  if [ $FUNGUS == 1 ]; then
    FUNGUS_PAR="--fungus"
  fi
  $THIS/../rm_all.sh data
  $THIS/../rm_all.sh info
  $THIS/../rm_all.sh output
  $THIS/../rm_all.sh run
  rm -f gmes.log
  rm -f run.cfg
  # Approximate
  /usr/local/gmes/4.69/bin/gmes_petap.pl  --ES  $FUNGUS_PAR  --sequence $ASM  --soft_mask 0  --cores $CORES  &>> $LOG
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
  $THIS/GeneMark2CDS $ASM genemark.gtf  -gtf  -gencode 1  -prot $ASM.prot  -prot_len_min 20  -ambig 10  -log $LOG  -qc
    # was: -prot_len_min 60

  section "$ASM.cds"
  $THIS/GeneMark2CDS $ASM genemark.gtf  -gtf  -gencode 1  -cds $ASM.cds  -prot_len_min $PROT_LEN  -complete  -log $LOG  -qc
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

if [ $UNIV ]; then
  section "prots2hmm_univ.sh"
  $THIS/prots2hmm_univ.sh $ASM $UNIV 0 $CORES $LOG >> $LOG
fi


gzip $ASM.prot


rm -f $LOG


