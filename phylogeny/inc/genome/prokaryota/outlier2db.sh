#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: Prok.id"
  echo "#2: Prok.outlier"
  exit 1
fi
ASM=$1
OUTLIER="$2"


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`

sqsh-ms -S $SERVER  -D $DATABASE  << EOT
  update Prok
    set outlier = '$OUTLIER'
    where id = $ASM
  if @@rowcount != 1
    raiserror ('update failed', 11, 1)
  go -m bcp 
EOT


