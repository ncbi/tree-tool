#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: Genome.id"
  echo "#2: GenomeTaxroot.taxroot_outlier"
  exit 1
fi


INC=$( dirname $0 )
SERVER=$( cat $INC/server )
DATABASE=$( cat $INC/database )
TAXROOT=$( cat $INC/../taxroot )

sqsh-ms -S $SERVER  -D $DATABASE  << EOT
  update GenomeTaxroot
    set taxroot_outlier = '$2'
    where     genome = $1
          and taxroot = $TAXROOT
  go -m bcp 
EOT
