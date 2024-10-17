#!/bin/bash --noprofile
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Compute dissimilarities for pairs of objects"
  echo "#1: input dissimilarity requests (pairs of objects) (absolute pathname)"
  echo "#2: new object file or directory, or ''. Object name is basename"
  echo "#3: output dissimilarities added to the pairs of objects (absolute pathname)"
  echo "#4: error log"
  exit 1
fi
REQUEST=$1
FILE_NEW="$2"
DISSIM=$3
LOG=$4


if [ $FILE_NEW ]; then
  error "FILE_NEW is not implemented"
fi


F=$THIS/query/$( uuidgen ).dissim
echo -e "$REQUEST\n$DISSIM\nEND" > $F
while [ -e $F ]; do
  sleep 0.001
done
while [ ! -s $DISSIM ]; do
  sleep 0.001
done


#rm -f $LOG
