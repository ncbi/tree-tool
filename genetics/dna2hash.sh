#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Invokes prodigal"
  echo "#1: prokaryotic genome DNA FASTA"
  echo "#2: output file with CDS hashes"
  echo "#3: file with target hashes | ''"
  echo "#4: output CDSs macthing target hashes | ''"
  exit 1
fi
FASTA=$1
OUT=$2
TARGET="$3"
CDS="$4"


if [ $TARGET ]; then
  $THIS/../check_file.sh $TARGET 1
fi


TMP=`mktemp`
comment $TMP

# F
GZ=`echo $FASTA | tr '.' '\n' | tail -1`
if [ $GZ == "gz" ]; then
  gunzip -c $FASTA > $TMP.fa
  F=$TMP.fa
else
  F=$FASTA
fi
prodigal  -o /dev/null  -i $F  -d $TMP.cds  &> /dev/null

$THIS/fasta2hash $TMP.cds $OUT  -cds  -gene_finder "prodigal"  -target_hashes $TARGET  -target_seq $TMP.ids > /dev/null

if [ $TARGET ]; then
 #$THIS/extractFastaDna $TMP.cds $TMP.ids -whole -qc > $CDS
  $THIS/filterFasta $TMP.cds  -target $TMP.ids  -whole -qc > $CDS
fi


#rm $TMP*

