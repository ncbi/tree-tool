#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../bash_common.sh
if [ $# -ne 8 ]; then
  echo "Create annotation files for #1 in the current directory"
  echo "Invokes GeneMark"
  echo "#1: eukaryotic genome name"
  echo "#2: eukaryotic DNA FASTA"
  echo "#3: Taxroot.id"
  echo "#4: universal HMM library (absolute path) or ''"
  echo "#5: Pfam HMM library (absolute path) or ''"
  echo "#6: use Pfam HMM cutoff (0/1)"
  echo "#7: number of cores"
  echo "#8: log file (absolute path)"
  echo "Time: #1 size: 10M - 1 hour; 4G - 6 days"
  exit 1
fi
ASM=$1
FASTA=$2
TAXROOT=$3
UNIV=$4
PFAM=$5
PFAM_CUTOFF=$6
CORES=$7
LOG=$8


#set -x


if [ ! -e $FASTA ]; then
  error "No $FASTA" >> $LOG
fi
if [ ! -s $FASTA ]; then
  error "Empty genome" >> $LOG
fi


function genemark_clean
{
  $THIS/../rm_all.sh data
  $THIS/../rm_all.sh info
  $THIS/../rm_all.sh output
  $THIS/../rm_all.sh run
  $THIS/../rm_all.sh prothint
  rm -f gmes.log run.cfg gmhmm.mod gmhmm_es.mod genemark_es.gtf
}


DIR=$( dirname $FASTA )

cd $DIR


ls -laF $FASTA
$THIS/dna2stat $FASTA  -log $LOG > $ASM.stat


if [ ! -e $ASM.prot ]; then
  if [ $TAXROOT == 58023 ]; then
    echo "dna2orfs 300" > annot_software
    $THIS/dna2prots $FASTA 1 300  -complexity_min 3.5  -no_x  -qc  -log $LOG > $ASM.prot
      # PAR  
    $THIS/fasta2hash $ASM.prot $ASM.hash-PRT -qc  -log $LOG
    ln -s $PWD/$ASM.hash-PRT $PWD/$ASM.hash-CDS
  else
    genemark_clean

    # http://topaz.gatech.edu/GeneMark/license_download.cgi: GeneMark-ES/ET/EP+ ver 4.73_lic
    #   fill in questionaire; press "I agree"; "Please donwload key: 64 bit"
    # $ gunzip gm_key_64.gz 
    # $ cp gm_key_64 ~/.gm_key 
    gmes_petap.pl | grep -w "version" 1> annot_software 2>> $LOG || true 

    FUNGUS_PAR=""
    if [ $TAXROOT == 4751 ]; then
      FUNGUS_PAR="--fungus"
        # https://pmc.ncbi.nlm.nih.gov/articles/PMC2593577/
    fi
    
    GMES_TYPE="--ES"
    if [ "$UNIV" ]; then
      UNIV_FA=$( echo $UNIV | sed 's/\.LIB$/.fa/1' )
      if [ -e $UNIV_FA ]; then
        GMES_TYPE="--EP  --dbep $UNIV_FA"
      fi
    fi
    echo "$GMES_TYPE" >> annot_software
    
    # genemark.gtf
    gmes_petap.pl  $GMES_TYPE  $FUNGUS_PAR  --sequence $FASTA  --soft_mask 0  --cores $CORES  &>> $LOG
    echo -e "\nAnnotation finished!\n" >> $LOG
    
    genemark_clean

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
fi


if [ "$PFAM" ]; then
  section "Pfam"
  $THIS/prots2hmm_hash.sh $ASM.prot $PFAM $PFAM_CUTOFF $ASM.HMM $ASM.hash-HMM $CORES $LOG >> $LOG
  gzip $ASM.HMM
fi

if [ "$UNIV" ]; then
  if [ "$TAXROOT" == 58023 ]; then  # PAR (Tracheophyta)
    section "tblastn2marker_euk.sh"
    UNIV_DIR=$( dirname $UNIV )
    $THIS/tblastn2marker_euk.sh $FASTA $UNIV_DIR/univ - $CORES $ASM.prot-univ $LOG
  else
    section "prots2hmm_univ.sh"
    $THIS/prots2hmm_univ.sh $ASM $UNIV 0 $CORES $LOG >> $LOG
  fi
fi


gzip $ASM.prot


rm -f $LOG


