#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
DB=GenomeHash
if [ $# -ne 2 ]; then
  echo "QC $DB.*"
  echo "Output: $DB.<SUF>.odd"
  echo "#1: genome/"
  echo "#2: #1 is large (0/1)"
  exit 1
fi
GENOME=$1
LARGE=$2


TMP=$( mktemp )


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


function qc_suf
{
  local SUF=$1
  if [ ! -e $DB.$SUF ]; then
    return
  fi
  section "$DB.$SUF"
  cut -f 1 $DB.$SUF | sort -u > $TMP.hash
  $THIS/../../setMinus $TMP $TMP.hash > $DB.$SUF.odd
  wc -l $DB.$SUF.odd
  if [ ! -s $DB.$SUF.odd ]; then
    rm $DB.$SUF.odd
  fi
}
qc_suf CDS
qc_suf PRT
qc_suf HMM


rm $TMP*
