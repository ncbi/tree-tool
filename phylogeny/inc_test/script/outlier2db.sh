#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Add an outlier to database"
  echo "#1: PDT.id"
  echo "#2: PDT.outlier"
  exit 1
fi

exit 0
