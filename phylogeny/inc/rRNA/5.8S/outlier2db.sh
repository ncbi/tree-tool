#!/bin/bash
source CPP_DIR/phylogeny/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: Locus5_8S.accession"
  echo "#2: Locus5_8S.outlier"
  exit 1
fi


INC=`dirname $0`
SERVER=`cat $INC/server`
DATABASE=`cat $INC/database`
sqsh-ms  -S $SERVER  -D $DATABASE  << EOT
  update Locus5_8S
    set outlier = '$2'
    where     accession = '$1'
          and selected = 1;
  go -m bcp 
EOT


