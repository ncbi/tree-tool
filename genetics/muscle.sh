#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "#1: input protein FASTA file"  
  echo "#2: output alignment"
  exit 1
fi
IN="$1"
OUT="$2"


muscle -sv -noanchors  -distance1 kbit20_3  -maxiters 2  -maxtrees 99  -in $IN  -out $OUT >& /dev/null
#clustalw -infile=$2 -outfile=$4.SEED -output=fasta -quicktree > /dev/null
