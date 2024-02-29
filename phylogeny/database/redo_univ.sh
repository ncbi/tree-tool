#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Input: hmm-univ.LIB, genome/#hash/#1/#1.prot.gz"
  echo "Output: genome/#hash/#1/#1.{univ,prot-univ}"
  echo "Update: log"
  echo "#1: assembly"
  exit 1
fi
ASM=$1


LOG=$PWD/log/$ASM

HASH=`$EXEC/file2hash $ASM`
cd genome/$HASH/$ASM/

if [ ! -e $ASM.prot ]; then
  gunzip $ASM.prot.gz
fi
$THIS/../genetics/prots2hmm_univ.sh $ASM ../../../hmm-univ.LIB 0 15 $LOG 
gzip $ASM.prot
