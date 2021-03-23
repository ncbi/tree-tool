#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: Genome.id"
  echo "#2: Genome.outlier"
  exit 1
fi
ASM=$1
OUTLIER="$2"


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`

sqsh-ms -S $SERVER  -D $DATABASE  << EOT
  update Genome
    set outlier = '$OUTLIER'
    where id = $ASM;
  go -m bcp 
EOT


