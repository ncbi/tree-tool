#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: Locus.id"
  echo "#2: Locus.outlier"
  exit 1
fi

sqsh-ms -S ""  -D uniColl  << EOT
  update Locus
    set outlier = '$2'
    where id = '$1';
  go -m bcp 
EOT


