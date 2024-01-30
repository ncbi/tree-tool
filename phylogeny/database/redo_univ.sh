#!/bin/bash --noprofile
THIS=`dirname $0`
EXEC=$THIS/../..
source $EXEC/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Output: genome/#hash/#1"
  echo "#1: assembly"
  exit 1
fi
ASM=$1


LOG=$PWD/log/$ASM

HASH=`$EXEC/file2hash $ASM`
cd genome/$HASH/$ASM/

gunzip $ASM.prot.gz
$EXEC/genetics/prots2hmm_univ.sh $ASM ../../../hmm.lib 0 $LOG 
gzip $ASM.prot
