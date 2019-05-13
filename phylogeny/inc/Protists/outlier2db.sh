#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: Genome.id"
  echo "#2: Genome.outlier"
  exit 1
fi

sqsh-ms -S PROTEUS  -D uniColl  << EOT
  update Genome
    set outlier = '$2'
    where id = $1;
  go -m bcp 
EOT


