#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Print assemblies with too few univesal proteins"
  echo "#1: assembly id. Input: genome/#1/#1.{prot-univ,hash}"
  echo "#2: Hashtype.id"
  echo "#3: Min. number of hashes"
  echo "#4: Max. number of hashes"
  exit 1
fi
ASM=$1
HASHTYPE=$2
MIN=$3
MAX=$4


H=`$THIS/../file2hash $ASM`
F=genome/$H/$ASM/$ASM.hash-$HASHTYPE

if [ ! -e $F ]; then
  error "$ASM No $F" 
fi

N=`cat $F | wc -l`
if [ $N -lt $MIN ]; then
  echo "$ASM small" 
  exit 0
fi
if [ $N -gt $MAX ]; then
  echo "$ASM large" 
  exit 0
fi



