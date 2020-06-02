#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Check the start/stop codons of a CDS FASTA file"
  echo "#1: CDS FASTA file"
  exit 1
fi
FASTA=$1


NAME=`basename $FASTA`

START_CODON=`tail -n +2 $FASTA | head -1 | cut -c 1-3`
SUFFIX=`echo $START_CODON | cut -c 2-3`
PREFIX=`echo $START_CODON | cut -c 1-2`
if [ "$PREFIX" != "at" -a "$SUFFIX" != "tg" ]; then
  echo $NAME
  exit 0
fi

STOP_CODON=`tail -n +2 $FASTA | tail -2 | tr '\n' ' ' | sed 's/ //g' | rev | cut -c 1-3 | rev`
if [ "$STOP_CODON" != "taa" -a "$STOP_CODON" != "tag" -a "$STOP_CODON" != "tga" ]; then
  echo $NAME
fi
