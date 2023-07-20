#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find closest genomes"
  echo "#1: Prok.id"
  echo "#2: directory or ''"
  exit 1
fi
GENOME=$1
DIR="$2"


INC=`dirname $0`

if [ -z $DIR ]; then
  DIR=$INC/../genome/$GENOME
fi


TMP=`mktemp`


blastp  -db $INC/prot-univ  -query $DIR/$GENOME.prot-univ  -comp_based_stats 0  -evalue 1e-10  -seg no  -max_target_seqs 100  -outfmt '6 sseqid' > $TMP
cut -f 1 -d '-' $TMP | sort | uniq -c | sort -n -r | awk '{print $2};' | grep -xv -f $INC/deleted.list | grep -xv $GENOME | head -100 || true


rm $TMP*
