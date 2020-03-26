#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: Locus.id"
  echo "#2: Locus.outlier"
  exit 1
fi

sqsh-ms -S PROTEUS  -D uniColl  << EOT
  update Locus
    set outlier = '$2'
    where id = '$1';
  if @@rowcount != 1
    raiserror ('Locus not found', 11, 1);
  go -m bcp 
EOT


