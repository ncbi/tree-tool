#!/bin/bash --noprofile
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find closest genomes among #3"
  echo "#1: object"
  echo "#2: new object directory or ''"
  echo "#3: subset of objects (absolute pathname)"
  echo "#4: output file (absolute pathname)"
  exit 1
fi
GENOME=$1
#DIR="$2"  # ??
SUBSET=$3
OUT=$4


F=$THIS/query/$( uuidgen ).closest
echo -e "$GENOME\n$SUBSET\n$OUT\nEND" > $F
while [ -e $F ]; do
  sleep 0.001
done


