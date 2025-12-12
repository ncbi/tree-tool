#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 3 ]; then
  echo "QC GenomeHash.*"
  echo "#1: genome/"
  echo "#2: #1 is large (0/1)"
  echo "#3: inc/"
  exit 1
fi
GENOME=$1
LARGE=$2
INC=$3


$THIS/../../check_file.sh $INC 0

SERVER=$( cat $INC/server )
DATABASE=$( cat $INC/database )
TAXROOT=$( cat $INC/../taxroot )


TMP=$( mktemp )
comment $TMP
#set -x


$THIS/../../check_file.sh $GENOME 0

# $TMP
if [ $LARGE -eq 1 ]; then
  $THIS/../../trav $GENOME "ls %d/%f" > $TMP
else
  ls $GENOME > $TMP
fi
wc -l $TMP
$THIS/../../sort.sh $TMP
sort -cu $TMP


sqsh-ms -S $SERVER  -D $DATABASE  -L bcp_rowsep=""  << EOT | sort > $TMP.annot
  select id
    from GenomeTaxroot_vw
    where     taxroot = $TAXROOT
          and prots is not null
  go -m bcp  
EOT
$THIS/../../setMinus $TMP.annot $TMP > $TMP.miss
wc -l $TMP.miss
if [ -s $TMP.miss ]; then
  error "Annotated genomes missing in $GENOME/"
fi


ls $GENOME.annot_failure.log/ | sed 's/\.gz$//1' > $TMP.log
$THIS/../../setIntersect.sh $TMP $TMP.log 0 > $TMP.genome-log
wc -l $TMP.genome-log
if [ -s $TMP.genome-log ]; then
  error "Genomes in $GENOME/ and $GENOME.annot_failure.log/"
fi


rm $TMP*
