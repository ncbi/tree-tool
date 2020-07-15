#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
UNKNOWN=0
NOVEL=999999
if [ $# != 3 ]; then
  echo "MLST typing"
  echo "Print: #1 {<ST-number>|$UNKNOWN|$NOVEL}, where $UNKNOWN means unknown, $NOVEL means novel"
  echo "#1: DNA identifier"
  echo "#2: DNA FASTA file: existing BLAST database is used or a temporary one is created"
  echo "#3: MLST directory with all.fa (all <locus>.fasta-files) and profiles.list"
  exit 1
fi
ID=$1
DNA=$2
MLST=$3


if [ ! $DNA ]; then
  echo "Empty DNA file name"
  exit 1
fi


TMP=`mktemp`  
#echo $TMP  


DB=$DNA
if [ ! -e $DB.nhr ]; then
  if [ ! -e $DNA ]; then
    echo "File $DNA does not exist"
    exit 1
  fi
  if [ ! -s $DNA ]; then
    echo "File $DNA is empty"
    exit 1
  fi
  DB=$TMP.db
  makeblastdb  -in $DNA  -out $DB  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
fi


ST=$UNKNOWN
N=`grep '^>' $MLST/all.fa | cut -f 1 -d '_' | sort -u | wc -l`
blastn  -db $DB  -query $MLST/all.fa  -dust no  -evalue 1e-50  -ungapped  -outfmt '6 qseqid qlen length nident' | awk '$2 == $3 && $2 == $4' | cut -f 1 | sort  > $TMP 
M1=`cut -f 1 -d '_' $TMP | sort | wc -l`
M2=`cut -f 1 -d '_' $TMP | sort -u | wc -l`
if [ $N -eq $M1 -a $M1 -eq $M2 ]; then
  L=`cut -f 2 -d '_' $TMP | tr '\n' '#' | sed 's/#$//1' | sed 's/^/#/1' | tr '#' '\t'`
  set +o errexit
  ST=`grep "$L$" $MLST/profiles.list | cut -f 1`
  S=$?
  set -o errexit
  if [ $S -eq 0 ]; then
    STS=`echo $ST | wc -l`
    if [ $STS -eq 1 ]; then
      if [ $ST -eq $UNKNOWN -o $ST -ge $NOVEL ]; then
        echo "$ST"
        exit 1
      fi
    fi
  else
    ST=$NOVEL
  fi
fi

echo "$ID $ST"


rm -f $TMP*  
