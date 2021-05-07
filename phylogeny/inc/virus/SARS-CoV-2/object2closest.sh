#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find closest virus"
  echo "#1: Virus.id"
  echo "#2: directory or ''"
  exit 1
fi
NEW_OBJ=$1
DIR="$2"


#set -o xtrace 


INC=`dirname $0`

if [ $DIR ]; then
  error "Novel object is not implemented"
fi


TMP=`mktemp`


SIZE=10000  # = ln(tree size)  * 1000
CPP_DIR/index_find.sh $INC/../mut.dna $INC/../mut.index $SIZE $INC/../mut.dna/$NEW_OBJ > $TMP
if [ -s $TMP ]; then
  grep -f $INC/../deleted.all -vx $TMP > $TMP.out
  head -100 $TMP.out
else
  if true; then  
    REAL_REF=`readlink $INC/../ref`
    basename $REAL_REF
  else
    echo "MW737421.1"
  fi
fi


rm $TMP*

