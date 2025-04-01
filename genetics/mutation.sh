#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../bash_common.sh
if [ $# -ne 10 ]; then
  echo "Print indels"
  echo "Cf. PUBMED 36677521"
  echo "#1: target DNA FASTA file"
  echo "#2: target name"
  echo "#3: reference DNA FASTA file"
  echo "#4: protein annotations in #2. Line format: <protein>\t<start>\t<stop> (1-based)"
  echo "#5: match_score, > 0"
	echo "#6: mismatch_score, > 0"
	echo "#7: gap_open, >= 0"
	echo "#8: gap_extent, > 0"
	echo "#9: indels only (0/1)"
  echo "#10: debug (0/1)"
  exit 1
fi
IN=$1
NAME=$2
REF=$3
ANNOT=$4
MATCH=$5
MISMATCH=$6
GAP_OPEN=$7
GAP_EXTENT=$8
INDEL=$9
DEBUG=${10}


TMP=$( mktemp )


ALIGN=""
if [ $DEBUG == 1 ]; then
  comment $TMP
  ALIGN="-alignment"
fi


$THIS/fasta2len $IN -qc -noprogress > $TMP.len
LEN=$( cut -f 2 $TMP.len )
if [ $LEN -lt 29600 -o $LEN -gt 31000 ]; then  # PAR
  echo -e "$NAME\tBAD: length $LEN"
  exit 
fi

$THIS/dna2stat $IN -noprogress > $TMP.stat
AMBIG=$( cut -f 3 $TMP.stat )
if (( $( echo "$AMBIG > 0.01" | bc -l ) )); then  # PAR
  echo -e "$NAME\tBAD: ambiguity $AMBIG"
  exit 
fi

# Unaligned target ends are not observed, probably because all targets have already been aligned to the reference before depositing to GenBank
$THIS/../dissim/seq2dissim $IN $REF  -mutation $TMP.mut  -ambig_start  -noambig 100  $ALIGN  -match_score $MATCH  -mismatch_score $MISMATCH  -gap_open $GAP_OPEN  -gap_extent $GAP_EXTENT > $TMP.dissim
if [ $INDEL == 1 ]; then
  grep '^INS' $TMP.mut >  $TMP.indel || true
  grep 'DEL$' $TMP.mut >> $TMP.indel || true
else
  mv $TMP.mut $TMP.indel
fi

$THIS/mutation_dna2prot $TMP.indel $REF $ANNOT  -annot_overlap > $TMP.out

sed 's/^/'$NAME'\t/1' $TMP.out

echo -e "$NAME\tDONE"


if [ $DEBUG == 0 ]; then
  rm $TMP*
fi

