#!/bin/bash --noprofile
EXEC=$PANFS/code
source $EXEC/cpp/bash_common.sh
ASM=$1


LOG=$PWD/log/$ASM

HASH=`$EXEC/cpp/file2hash $ASM`
cd genome/$HASH/$ASM/

gunzip $ASM.prot.gz
$EXEC/cpp/genetics/prots2hmm_univ.sh $ASM ../../../hmm.lib 0 $LOG 
gzip $ASM.prot
