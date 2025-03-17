#!/bin/bash --noprofile
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Record an outlier in a database"
  echo "#1: object id"
  echo "#2: outlier type"
  exit 1
fi
OBJ=$1
OUTLIER=$2


#echo -e "${OBJ}\t${OUTLIER}" >> $THIS/outlier
exit 1

