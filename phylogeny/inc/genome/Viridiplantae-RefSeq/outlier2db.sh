#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: RGenome.acc_ver"
  echo "#2: RGenome.outlier"
  exit 1
fi
ACC_VER=$1
OUTLIER=$2


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`

sqsh-ms -S $SERVER  -D $DATABASE  << EOT
  update RGenome
    set outlier = '$OUTLIER'
    where acc_ver = '$ACC_VER'
  go -m bcp 
EOT
