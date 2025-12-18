#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: Genome.id"
  echo "#2: GenomeTaxroot.outlier"
  exit 1
fi
GENOME=$1
OUTLIER="$2"


INC=$( dirname $0 )
SERVER=$( cat $INC/server )
DATABASE=$( cat $INC/database )
TAXROOT=$( cat $INC/../taxroot )

sqsh-ms -S $SERVER  -D $DATABASE  << EOT
  update GenomeTaxroot
    set outlier = '$OUTLIER'
    where     genome = $GENOME
          and taxroot = $TAXROOT
  go -m bcp 
EOT
