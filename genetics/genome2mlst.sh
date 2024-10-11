#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
# Cf. blastn2mlst.cpp
UNKNOWN=0
NOVEL=999999
#
if [ $# != 4 ]; then
  echo "MLST typing"
  echo "Print: #1 {<ST-number>|$UNKNOWN|$NOVEL}, where $UNKNOWN means unknown, $NOVEL means novel"
  echo "#1: DNA identifier"
  echo "#2: DNA FASTA file: existing BLAST database is used or a temporary one is created"
  echo "#3: MLST directory with all.fa (all <locus>.fasta-files) and profiles.list"
  echo "#3: verbose (0/1)"
  exit 1
fi
ID=$1
DNA=$2
MLST=$3
VERB=$4


$THIS/../check_file.sh "$DNA" 1

if [ ! -s $DNA ]; then
  error "File $DNA is empty"
fi


TMP=$( mktemp )
#comment $TMP  


DB=$DNA
if [ ! -e $DB.nhr ]; then
  DB=$TMP.db
  makeblastdb  -in $DNA  -out $DB  -dbtype nucl  -blastdb_version 4  -logfile /dev/null 
fi

# $TMP.blastn
blastn  -db $DB  -query $MLST/all.fa  -dust no  -evalue 1e-50  -ungapped  -outfmt '6 qseqid sseqid sstart send qlen length nident' | awk '$5 == $6 && $5 == $7' | sort > $TMP.blastn
#                                                                                    1      2      3      4    5    6      7
if [ $VERB == 1 ]; then
  echo -e "qseqid\tsseqid\tsstart\tsend\tqlen\tlength\tnident"
  cat $TMP.blastn
fi

ST=$( $THIS/blastn2mlst $TMP.blastn $MLST )
echo -e "$ID\t$ST"


rm -f $TMP*  
