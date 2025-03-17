#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find closest objects among #3"
  echo "#1: object id"
  echo "#2: new object directory | ''"
  echo "#3: subset of objects id's (absolute pathname)"
  echo "#4: output file (absolute pathname)"
  exit 1
fi
OBJ=$1
DIR=$2
SUBSET=$3
OUT=$4


INC=$( dirname $0 )

if [ -z $DIR ]; then
  DIR=$INC/../seq
fi

CPP_DIR/genetics/prot2closest.sh $DIR/$OBJ $INC/../db $SUBSET > $OUT

