#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: Locus.id"
  echo "#2: Locus.outlier"
  exit 1
fi


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
sqsh-ms  -S $SERVER  -D $DATABASE  << EOT
  update Locus
    set outlier = '$2'
    where id = '$1';
  go -m bcp 
EOT